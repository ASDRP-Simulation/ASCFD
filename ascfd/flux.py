
from onedim.euler import Euler
from onedim.reconstruct import weno5_reconstruction
from onedim.constants import *
from onedim.grid import Grid1D
from onedim.constants import Constants

import numpy as np
import sys

class Flux:

    def __init__(self, a_constants: Constants, a_type: str):

        self.type = a_type
        self.c = a_constants

        self.euler = Euler(self.c)


        if self.type == "LF":
            self.flux_method = self.lax_friedrichs
        else:
            raise RuntimeError(f"Flux method not supported: {self.type}")


    def getFlux(self, a_grid):
        return self.flux_method(a_grid)


    def lax_friedrichs(self, a_grid):
        a_grid.assert_variable_type("prim")

        #get density 
        if self.c.system == "euler2D":
            density = a_grid.grid[self.c.RHOCOMP]
        else:
            raise RuntimeError("Density method needs to be implemented.")

                
        a = np.sqrt(gamma * a_grid.grid[PCOMP] / density)

        max_speed = np.max(np.abs(a_grid.grid[self.c.UCOMP]) + a)



        #eq 29
        sum = 0
        for i in range(self.c.NS):
            sum += (a_grid.grid[self.c.RHOCOMP + i] / self.c.mW[i])
        mBar = density / sum

        #eq 4
        sum = 0
        Y_is = []
        for i in range(self.c.NS):

            #Y = rhoY/rho
            Y_i = a_grid.grid[self.c.RHOCOMP + i]/density
            Y_is.append(Y_i)

            sum += (1 / (1 - self.c.gammas[i])) * (Y_i / self.c.mW[i])

        G = mBar * sum
        gammaBar = 1/G + 1



        #don't reconstruct cell faces just use FD methods
        U = a_grid.grid

        # conservative variables 
        consU = self.euler.prim_to_cons(U)

        # analytical flux
        f = self.euler.flux(U)

        #pointwise LF Flux
        numFluxL = np.zeros_like(a_grid.grid)
        numFluxR = np.zeros_like(a_grid.grid)

        # Computes LF Flux
        # Loop through all the cells except for the outermost ghost cells.
        for i in range(1, len(a_grid.x) - 1):

            leftWave = np.abs(a_grid.grid[self.c.UCOMP, i]) + np.sqrt(gammaBar[i] * a_grid.grid[self.c.PCOMP,i] / density[i])
            rightWave = np.abs(a_grid.grid[self.c.UCOMP, i+1]) + np.sqrt(gammaBar[i+1] * a_grid.grid[self.c.PCOMP,i+1] / density[i+1])
            sMaxR = max(leftWave, rightWave)


            leftWave = np.abs(a_grid.grid[self.c.UCOMP, i-1]) + np.sqrt(gammaBar[i-1] * a_grid.grid[self.c.PCOMP,i-1] / density[i-1])
            rightWave = np.abs(a_grid.grid[self.c.UCOMP, i]) + np.sqrt(gammaBar[i] * a_grid.grid[self.c.PCOMP,i] / density[i])
            sMaxL = max(leftWave, rightWave)

            for icomp in range(self.c.NUMQ):

                # Compute the Lax-Friedrichs flux                 
                numFluxR[icomp, i] = 0.5*(f[icomp, i+1] + f[icomp, i]) - 0.5*sMaxR * (consU[icomp,i+1] - consU[icomp, i])
                numFluxL[icomp, i] = 0.5*(f[icomp, i] + f[icomp,i-1]) - 0.5*sMaxL * (consU[icomp,i] - consU[icomp, i-1])

        return consU, numFluxR, numFluxL

