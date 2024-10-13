from onedim.grid import Grid1D
from onedim.constants import Constants
import onedim.ics as ics
import glob
import matplotlib.animation as animation
import matplotlib.ticker as ticker
import copy

import sys
from onedim.flux import Flux


from onedim.euler import Euler

from onedim.reconstruct import weno5_reconstruction
from onedim.bcs import BoundaryConditions

import numpy as np
import matplotlib.pyplot as plt
import os

class Simulation:
    def __init__(self, a_inputs):
        self.inp = a_inputs
        self.c = Constants(a_inputs)
        self.euler = Euler(self.c)

        self.grid = Grid1D(self.inp.xlim, self.inp.nx, self.inp.numghosts, self.c.NUMQ)
        self.bcs = BoundaryConditions(self.grid, self.inp.bc_lo, self.inp.bc_hi)
        self.flux = Flux(self.c, self.inp.flux)

        self.applyICS()

        self.bcs.apply_bcs()
        self.grid.check_grid(self.c)

        self.t = self.inp.t0
        self.timestepNum = 0

 
        #-1 is no output. Always output ICs if we are outputting.
        if self.inp.output_freq >= 0:
            self.output()

    def run(self):
        # os.makedirs('simulation_frames', exist_ok=True)


        while (self.t < self.inp.t_finish) and self.timestepNum < self.inp.nt:
            print("Timestep: ", self.timestepNum, "  Current time: ", self.t)

            self.bcs.apply_bcs()


            self.grid.assert_variable_type("prim")

            
            #get density 
            if self.inp.system == "euler1D":
                density = self.grid.grid[self.c.RHOCOMP]
            
            elif self.inp.system == "euler1DNS2" or self.inp.system == "4HRM":
                density = np.zeros_like(self.grid.grid[0])
                
                #rho = rhoY1 + rhoY2
                for i in range(self.c.NS):
                    density += self.grid.grid[i]
            else:
                raise RuntimeError("Density method needs to be implemented.")


            a = np.sqrt(self.c.gamma * self.grid.grid[self.c.PCOMP] / density)
            max_speed = np.max(np.abs(self.grid.grid[self.c.UCOMP]) + a)
            dt = min(self.inp.cfl * self.grid.dx / max_speed, self.inp.t_finish - self.t)
            
            U_new = np.zeros_like(self.grid.grid) #/ 0  # np.nans_like lol
            
            
            if  self.inp.timeStepper == "RK1":
                #returns numerical flux and conservative varaibles at interface
                consP, numericalFluxP, numericalFluxM = self.flux.getFlux(self.grid)

                for i in range(self.grid.Nghost, self.grid.Nx + self.grid.Nghost):
                    for icomp in range(self.c.NUMQ):
                        U_new[icomp, i] = consP[icomp, i] - (dt / self.grid.dx) * (
                            numericalFluxP[icomp, i] - numericalFluxM[icomp, i]
                        )


            elif self.inp.timeStepper == "RK4":

                rhs = np.zeros_like(self.grid.grid) #/ 0  # np.nans_like lol
                self.grid.check_grid(self.c,prim=True)

                u_start = self.euler.prim_to_cons(self.grid.grid)


                consP, numericalFluxP, numericalFluxM = self.flux.getFlux(self.grid)
                for i in range(self.grid.Nghost, self.grid.Nx + self.grid.Nghost):
                    for icomp in range(self.c.NUMQ):
                        rhs[icomp, i] =  (1 / self.grid.dx) * (
                            numericalFluxP[icomp, i] - numericalFluxM[icomp, i]
                        )

                k1 = dt*rhs
                self.grid.grid = u_start - k1 / 2

                self.bcs.apply_bcs()
                self.grid.check_grid(self.c,cons=True)


                self.grid.grid = self.euler.cons_to_prim(u_start - k1 / 2)
                self.bcs.apply_bcs()
                self.grid.check_grid(self.c,prim=True)


                consP, numericalFluxP, numericalFluxM = self.flux.getFlux(self.grid)
                for i in range(self.grid.Nghost, self.grid.Nx + self.grid.Nghost):
                    for icomp in range(self.c.NUMQ):
                        rhs[icomp, i] =  (1 / self.grid.dx) * (
                            numericalFluxP[icomp, i] - numericalFluxM[icomp, i]
                        )
                
                k2 = dt * rhs
                self.grid.grid = self.euler.cons_to_prim(u_start - k2 / 2)
                self.bcs.apply_bcs()
                self.grid.check_grid(self.c,prim=True)


                consP, numericalFluxP, numericalFluxM = self.flux.getFlux(self.grid)
                for i in range(self.grid.Nghost, self.grid.Nx + self.grid.Nghost):
                    for icomp in range(self.c.NUMQ):
                        rhs[icomp, i] =  (1 / self.grid.dx) * (
                            numericalFluxP[icomp, i] - numericalFluxM[icomp, i]
                        )

                k3 = dt*rhs
                self.grid.grid = self.euler.cons_to_prim(u_start - k3)
                self.bcs.apply_bcs()
                #self.grid.check_grid(self.c)





                consP, numericalFluxP, numericalFluxM = self.flux.getFlux(self.grid)
                for i in range(self.grid.Nghost, self.grid.Nx + self.grid.Nghost):
                    for icomp in range(self.c.NUMQ):
                        rhs[icomp, i] =  (1 / self.grid.dx) * (
                            numericalFluxP[icomp, i] - numericalFluxM[icomp, i]
                        )
                k4 = dt*rhs
                U_new = u_start - ((1/6)*k1  + (1/3)*k2 + (1/3)*k3 + (1/6)*k4)
                self.bcs.apply_bcs()
                self.grid.check_grid(self.c)


            else:
                raise RuntimeError("Timestepping method not supported.")



            self.grid.set(U_new)
            self.grid.transform(self.euler.cons_to_prim, "prim")
            
            self.bcs.apply_bcs()


            self.timestepNum += 1
            self.t += dt

            #always output the last timestep.
            if (self.timestepNum % self.inp.output_freq == 0) or (self.timestepNum == self.inp.nt-1):
                
                self.output()


            # DEBUG
            self.grid.check_grid(self.c)

        if self.inp.make_movie:
            self.generate_movie()
    
        print("SUCCESS!")
        return self.grid

    def plot(self):
        if not os.path.exists(self.inp.output_dir):
            os.makedirs(self.inp.output_dir)

        fig, axs = plt.subplots(3, 1, figsize=(10, 15))

        axs[0].scatter(self.grid.x, self.grid.grid[self.c.RHOCOMP, :], c="black")
        axs[0].set_ylabel("Density")

        axs[1].scatter(self.grid.x, self.grid.grid[self.c.UCOMP, :],  c="black")
        axs[1].set_ylabel("Velocity")

        axs[2].scatter(self.grid.x, self.grid.grid[self.c.PCOMP, :],  c="black")
        axs[2].set_ylabel("Pressure")

        axs[0].set_title(f"Time: {self.t:.4f}")
        plt.savefig(f"{self.inp.output_dir}/plot_dt{str(self.timestepNum).zfill(6)}")
        plt.close()

    

    def applyICS(self):

        if self.inp.system == "euler1D":
            if self.inp.ics == "sodshocktube":
                self.grid.fill_grid(ics.sod_shock_tube)
            elif self.inp.ics == "shuoshershocktube":
                self.grid.fill_grid(ics.shu_osher_shock_tube)
            elif self.inp.ics == "lax_problem":
                self.grid.fill_grid(ics.lax_problem)
            else:
                raise RuntimeError("[FLUID] ICS not valid.")
            
        elif self.inp.system  == "euler1DNS2":
            if self.inp.ics == "NS2_smooth_interfaces":
                self.grid.fill_grid(ics.NS2_smooth_interfaces)
        else:
            raise RuntimeError("[FLUID] ICS not valid.")

        

    def output(self):
        if not os.path.exists(self.inp.output_dir):
            os.makedirs(self.inp.output_dir)
        
        # File naming convention: output_timestepNum.txt
        output_filename = os.path.join(self.inp.output_dir, f"output_{str(self.timestepNum).zfill(6)}.txt")
        output_plotname = os.path.join(self.inp.output_dir, f"output_{str(self.timestepNum).zfill(6)}.png")

        with open(output_filename, 'w') as f:
            # Write header
            f.write(f"# Time: {self.t:.4f}\n")
            f.write("# x, density, velocity, pressure\n")
            
            for i in range(len(self.grid.x)):
                x = self.grid.x[i]
                components = [self.grid.grid[q, i] for q in range(self.c.NUMQ)]
                f.write(f"{x:.12f}, " + ", ".join(f"{comp:.8f}" for comp in components) + "\n")
        
        if self.c.NUMQ == 3:
            fig, axs = plt.subplots(3, 1, figsize=(10, 15))
        elif self.c.NUMQ == 4:
            fig, axs = plt.subplots(2,2, figsize=(15, 15))
        else:
            print("self.c.NUMQ: ", self.c.NUMQ)
            raise RuntimeError("System not implemented for movie.")
        
        axs = axs.ravel()  # Flatten the array to index by i

    
        for i in range(self.c.NUMQ):
            #axs[i].clear()
            axs[i].scatter(self.grid.x, self.grid.grid[i], c="black")
            axs[i].set_ylabel(self.c.variable_names[i])


        axs[0].set_title(f"Time: {self.t:.4f}, Timestep: {self.timestepNum}")

        # for ax in axs.flat:
        #     ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useOffset=False))

        # axs[2].set_ylim((0.99, 1.01))
        # axs[3].set_ylim((0.89, .91))

        fig.savefig(output_plotname)
        plt.close()

    





    def generate_movie(self):
        # Create a directory for the frames if it doesn't exist
        frames_dir = os.path.join(self.inp.output_dir, "frames")
        if not os.path.exists(frames_dir):
            os.makedirs(frames_dir)

        # List all the output files and sort them
        #output_files = sorted(glob.glob(os.path.join(self.inp.output_dir, "output_*.png")))
            
        #movie_filename = os.path.join(self.inp.output_dir, "simulation_movie.mp4")


        ffmpeg_command = f"ffmpeg -y -framerate 24 -i {self.inp.output_dir}/output_%06d.png -c:v libx264 -pix_fmt yuv420p {self.inp.output_dir}/movie.mp4"
        os.system(ffmpeg_command)



        # if self.c.NUMQ == 3:
        #     fig, axs = plt.subplots(3, 1, figsize=(10, 15))
        # elif self.c.NUMQ == 4:
        #     fig, axs = plt.subplots(2,2, figsize=(10, 15))
        # else:
        #     print("self.c.NUMQ: ", self.c.NUMQ)
        #     raise RuntimeError("System not implemented for movie.")
        
        # axs = axs.ravel()  # Flatten the array to index by i




        # def update_plot(file):
        #     data = np.loadtxt(file, delimiter=',', skiprows=2)
        #     x = data[:, 0]

        #     #data is indexed by
        #     # data[:,i] where i = 0 for x, i = 1 for icomp1, i=2 for icomp2

        #     with open(file, 'r') as f:
        #         lines = f.readlines()
        #         time_line = lines[0]
        #         time = float(time_line.split(':')[1].strip())

        #     timestep = int(file.split('_')[-1].split('.')[0])


        #     for i in range(self.c.NUMQ):
        #         axs[i].clear()
                

        #         axs[i].scatter(x, data[:,i+1], c="black")
        #         axs[i].set_ylabel(self.c.variable_names[i])


        #     axs[0].set_title(f"Time: {time:.4f}, Timestep: {timestep}")
        
        # # Create an animation by updating the plot for each output file
        # ani = animation.FuncAnimation(fig, update_plot, frames=output_files, repeat=False)

        # # Save the animation as a movie file using ffmpeg
        # movie_filename = os.path.join(self.inp.output_dir, "simulation_movie.mp4")
        # ani.save(movie_filename, writer='ffmpeg', fps=10)

        # print(f"Movie saved as {movie_filename}")
