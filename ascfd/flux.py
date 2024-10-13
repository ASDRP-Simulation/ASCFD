
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


        if self.type == "weno5_js":
            self.flux_method = self.weno_js
        elif self.type == "simple_flux":
            self.flux_method = self.simple_flux
        elif self.type == "PE_LLF":
            #pressure equilbrium local lax fried.
            self.flux_method = self.PE_MP_flux
        else:
            raise RuntimeError(f"Flux method not supported: {self.type}")


    def getFlux(self, a_grid):
        return self.flux_method(a_grid)

    # def weno_js(self, a_grid):
    #     a_grid.assert_variable_type("prim")

    #     #get density 
    #     if self.c.system == "euler1D":
    #         density = a_grid.grid[self.c.RHOCOMP]
        
    #     elif self.c.system == "euler1DNS2":
    #         density = np.zeros_like(a_grid.grid[0])
            
    #         #rho = rhoY1 + rhoY2
    #         for i in range(self.c.NS):
    #             density += a_grid.grid[i]
    #     else:
    #         raise RuntimeError("Density method needs to be implemented.")

    #     a = np.sqrt(gamma * a_grid.grid[PCOMP] / density)
    #     max_speed = np.max(np.abs(a_grid.grid[self.c.UCOMP]) + a)

    #     # reconstruct primitive at i+/-1/2 cell face
    #     primP, primM = weno5_reconstruction(a_grid.grid, a_grid)

    #     # conservative variables at cell interface.
    #     consP = self.euler.prim_to_cons(primP)
    #     consM = self.euler.prim_to_cons(primM)

    
    #     # compute analytical flux at cell interface.
    #     fR = self.euler.flux(primP) 
    #     fL = self.euler.flux(consM)


    #     #i think we actually want to use the cell center values for the second half of the 
    #     #LF flux.
    #     cellCenterCons = self.euler.prim_to_cons(a_grid.grid)


    #     #P/R [i] = i+1/2
    #     #M/L [i] = i-1/2


 
    #     LFFluxP = np.zeros_like(a_grid.grid)
    #     LFFluxM = np.zeros_like(a_grid.grid)

    #     # Computes LF Flux
    #     # Loop through all the cells except for the outermost ghost cells.
    #     for i in range(1, len(a_grid.x) - 1):
    #         for icomp in range(self.c.NUMQ):
    #             # Compute the Lax-Friedrichs flux at i+1/2 and i-1/2 interfaces
    #             # LFFluxP[icomp, i] = 0.5 * (
    #             #     fR[icomp, i ] + fR[icomp, i-1]
    #             # ) - 0.5 * max_speed * (cellCenterCons[icomp, i+1 ] - cellCenterCons[icomp, i])
                
                
    #             # LFFluxM[icomp, i] = 0.5 * (
    #             #     fL[icomp, i ] + fL[icomp, i+1]
    #             # ) - 0.5 * max_speed * (cellCenterCons[icomp, i ] - cellCenterCons[icomp, i-1])


    #             LFFluxP[icomp,i] = fR[icomp,i]
    #             LFFluxM[icomp,i] = fL[icomp,i]


    #     return consP, LFFluxP, LFFluxM

    def weno_js(self, a_grid):
        a_grid.assert_variable_type("prim")

        a = np.sqrt(gamma * a_grid.grid[PCOMP] / a_grid.grid[RHOCOMP])
        max_speed = np.max(np.abs(a_grid.grid[UCOMP]) + a)

        # reconstruct primitive at i+1/2 cell face
        primPR, primPL = weno5_reconstruction(a_grid.grid, a_grid)
        # primP = SEDAS_apriori(grid.grid, grid)

        # conservative variables at cell interface.
        # print("shape:" , np.shape(primP))
        consP = self.euler.prim_to_cons(primPR)
        consPL = self.euler.prim_to_cons(primPL)

        # compute analytical flux at cell interface.
        fR = self.euler.flux(primPR)

        U_new = np.ones_like(a_grid.grid) #/ 0  # np.nans_like lol
        LFFluxP = np.zeros_like(a_grid.grid)
        LFFluxM = np.zeros_like(a_grid.grid)

        # Computes LF Flux
        # Loop through all the cells except for the outermost ghost cells.
        for i in range(1, len(a_grid.x) - 1):
            for icomp in range(NUMQ):
                # Compute the Lax-Friedrichs flux at i+1/2 and i-1/2 interfaces
                LFFluxP[icomp, i] = 0.5 * (
                    fR[icomp, i + 1] + fR[icomp, i]
                ) - 0.5 * max_speed * (consP[icomp, i + 1] - consP[icomp, i])

                LFFluxM[icomp, i] = 0.5 * (
                    fR[icomp, i ] + fR[icomp, i-1]
                ) - 0.5 * max_speed * (consP[icomp, i ] - consP[icomp, i-1])

        return consP, LFFluxP, LFFluxM


    def simple_flux(self, a_grid):
        a_grid.assert_variable_type("prim")

        #get density 
        if self.c.system == "euler1D":
            density = a_grid.grid[self.c.RHOCOMP]

            for i in range(len(density)):
                if density[i] < 0:
                    print("negative density")
                    sys.exit()
        
        elif self.c.system == "euler1DNS2":
            density = np.zeros_like(a_grid.grid[0])
            
            #rho = rhoY1 + rhoY2
            for i in range(self.c.NS):
                density += a_grid.grid[i]
        else:
            raise RuntimeError("Density method needs to be implemented.")

                
        a = np.sqrt(gamma * a_grid.grid[PCOMP] / density)
        # for i in range(len(a)):

        #     if gamma * a_grid.grid[PCOMP,i] / density[i] <= 0:
        #         print("bad")
        #         print("gamma: ", gamma)
        #         print(a_grid.grid[PCOMP,i])
        #         print(density[i])
        #         print("i; ", i)
        #         sys.exit()

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


    def PE_MP_flux(self, a_grid):
        '''
        Fully conservative and pressure-equilibrium preserving
        scheme for compressible multi-component flows.
        No shock capturing!!!
        '''
         
        a_grid.assert_variable_type("prim")

        if self.c.system != "euler1DNS2":
            raise RuntimeError("PE_MP_flux flux is only for euler1DNS2.")

        density = np.zeros_like(a_grid.grid[0])
        #rho = rhoY1 + rhoY2
        for i in range(self.c.NS):
            density += a_grid.grid[self.c.RHOCOMP + i]

       
        # a = np.sqrt(gamma * a_grid.grid[PCOMP] / density)
        # max_speed = np.max(np.abs(a_grid.grid[self.c.UCOMP]) + a)



        #eq 29
        mBar = density / ( a_grid.grid[self.c.RHOCOMP] / self.c.mW[0] + a_grid.grid[self.c.RHOCOMP + 1] / self.c.mW[1])

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

        # conservative variables 
        consU = self.euler.prim_to_cons(a_grid.grid)


        #pointwise LF Flux
        numFluxL = np.zeros_like(a_grid.grid)
        numFluxR = np.zeros_like(a_grid.grid)

        # Computes LF Flux
        # Loop through all the cells except for the outermost ghost cells.
        
        C1 = np.zeros_like(a_grid.grid[0])
        C2 = np.zeros_like(a_grid.grid[0])
        Mu = np.zeros_like(a_grid.grid[0])
        pi = np.zeros_like(a_grid.grid[0])
        K = np.zeros_like(a_grid.grid[0])
        I = np.zeros_like(a_grid.grid[0])
        P = np.zeros_like(a_grid.grid[0])

        for i in range(1, len(a_grid.x) - 1):
            #for icomp in range(self.c.NUMQ):

                #just computing F_i + 1/2
            
            #eq 35
            phi_plus = (mBar[i+1]/density[i+1]) * (density[i]/mBar[i])
            phi_minus = (mBar[i]/density[i]) * (density[i+1]/mBar[i+1])
            
            
            #table 1
            C1[i] = (phi_minus * a_grid.grid[self.c.RHOCOMP, i] + phi_plus * a_grid.grid[self.c.RHOCOMP, i+1])/2  * (a_grid.grid[self.c.UCOMP, i] + a_grid.grid[self.c.UCOMP, i+1])/2
            C2[i] = (phi_minus * a_grid.grid[self.c.RHOCOMP+1, i] + phi_plus * a_grid.grid[self.c.RHOCOMP+1, i+1])/2  * (a_grid.grid[self.c.UCOMP, i] + a_grid.grid[self.c.UCOMP, i+1])/2

            Mu[i] = (phi_minus * density[i] + phi_plus * density[i+1])/2 \
                  * (a_grid.grid[self.c.UCOMP, i] + a_grid.grid[self.c.UCOMP, i+1])/2 \
                  * (a_grid.grid[self.c.UCOMP, i] + a_grid.grid[self.c.UCOMP, i+1])/2

            pi[i] = (a_grid.grid[self.c.PCOMP, i] + a_grid.grid[self.c.PCOMP, i+1])/2

            K[i] = (phi_minus * density[i] + phi_plus*density[i+1])/2 \
                 * (a_grid.grid[self.c.UCOMP, i] + a_grid.grid[self.c.UCOMP, i+1])/2 \
                 * (a_grid.grid[self.c.UCOMP, i] * a_grid.grid[self.c.UCOMP, i+1])/2

            I[i] = (a_grid.grid[self.c.PCOMP, i] * G[i] + a_grid.grid[self.c.PCOMP, i+1] * G[i+1])/2 \
                 * (a_grid.grid[self.c.UCOMP,i] + a_grid.grid[self.c.UCOMP,i+1])/2

            P[i] = (a_grid.grid[self.c.UCOMP,i] * a_grid.grid[self.c.PCOMP, i+1]  + a_grid.grid[self.c.UCOMP,i+1] * a_grid.grid[self.c.PCOMP, i])/2

        for i in range(1, len(a_grid.x) - 1):

            numFluxR[self.c.RHOCOMP, i] = C1[i]
            numFluxL[self.c.RHOCOMP, i] = C1[i-1]

            numFluxR[self.c.RHOCOMP + 1, i] = C2[i]
            numFluxL[self.c.RHOCOMP + 1, i] = C2[i-1]

            numFluxR[self.c.UCOMP, i] = Mu[i] + pi[i]
            numFluxL[self.c.UCOMP, i] = Mu[i-1] + pi[i-1]

            numFluxR[self.c.PCOMP, i] =  K[i] + I[i] + P[i]
            numFluxL[self.c.PCOMP, i] =  K[i-1] + I[i-1] + P[i-1]

        return consU, numFluxR, numFluxL



    def PE_MP_flux_LF_AD(self, a_grid):
        '''
        Fully conservative and pressure-equilibrium preserving
        scheme for compressible multi-component flows.

        Combined with the ideas in:
        Stable, entropy-consistent, and localized artificial-viscosity method for capturing shocks and contact discontinuities

        To use the pressure equilbirum fluxes as the non dissipative part in LF but have some dissipation.
        '''
         
        a_grid.assert_variable_type("prim")

        if self.c.system != "euler1DNS2":
            raise RuntimeError("PE_MP_flux flux is only for euler1DNS2.")

        density = np.zeros_like(a_grid.grid[0])
        #rho = rhoY1 + rhoY2
        for i in range(self.c.NS):
            density += a_grid.grid[self.c.RHOCOMP + i]


        a = np.sqrt(gamma * a_grid.grid[PCOMP] / density)
        max_speed = np.max(np.abs(a_grid.grid[self.c.UCOMP]) + a)

       
        # a = np.sqrt(gamma * a_grid.grid[PCOMP] / density)
        # max_speed = np.max(np.abs(a_grid.grid[self.c.UCOMP]) + a)



        #eq 29
        mBar = density / ( a_grid.grid[self.c.RHOCOMP] / self.c.mW[0] + a_grid.grid[self.c.RHOCOMP + 1] / self.c.mW[1])

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

        # conservative variables 
        consU = self.euler.prim_to_cons(a_grid.grid)


        #pointwise LF Flux
        numFluxL = np.zeros_like(a_grid.grid)
        numFluxR = np.zeros_like(a_grid.grid)

        pressureEquilbriumFlux = np.zeros_like(a_grid.grid)


        # Computes LF Flux
        # Loop through all the cells except for the outermost ghost cells.
        
        C1 = np.zeros_like(a_grid.grid[0])
        C2 = np.zeros_like(a_grid.grid[0])
        Mu = np.zeros_like(a_grid.grid[0])
        pi = np.zeros_like(a_grid.grid[0])
        K = np.zeros_like(a_grid.grid[0])
        I = np.zeros_like(a_grid.grid[0])
        P = np.zeros_like(a_grid.grid[0])

        for i in range(1, len(a_grid.x) - 1):
            #for icomp in range(self.c.NUMQ):

                #just computing F_i + 1/2
            
            #eq 35
            phi_plus = (mBar[i+1]/density[i+1]) * (density[i]/mBar[i])
            phi_minus = (mBar[i]/density[i]) * (density[i+1]/mBar[i+1])
            
            
            #table 1
            C1[i] = (phi_minus * a_grid.grid[self.c.RHOCOMP, i] + phi_plus * a_grid.grid[self.c.RHOCOMP, i+1])/2  * (a_grid.grid[self.c.UCOMP, i] + a_grid.grid[self.c.UCOMP, i+1])/2
            C2[i] = (phi_minus * a_grid.grid[self.c.RHOCOMP+1, i] + phi_plus * a_grid.grid[self.c.RHOCOMP+1, i+1])/2  * (a_grid.grid[self.c.UCOMP, i] + a_grid.grid[self.c.UCOMP, i+1])/2

            Mu[i] = (phi_minus * density[i] + phi_plus * density[i+1])/2 \
                  * (a_grid.grid[self.c.UCOMP, i] + a_grid.grid[self.c.UCOMP, i+1])/2 \
                  * (a_grid.grid[self.c.UCOMP, i] + a_grid.grid[self.c.UCOMP, i+1])/2

            pi[i] = (a_grid.grid[self.c.PCOMP, i] + a_grid.grid[self.c.PCOMP, i+1])/2

            K[i] = (phi_minus * density[i] + phi_plus*density[i+1])/2 \
                 * (a_grid.grid[self.c.UCOMP, i] + a_grid.grid[self.c.UCOMP, i+1])/2 \
                 * (a_grid.grid[self.c.UCOMP, i] * a_grid.grid[self.c.UCOMP, i+1])/2

            I[i] = (a_grid.grid[self.c.PCOMP, i] * G[i] + a_grid.grid[self.c.PCOMP, i+1] * G[i+1])/2 \
                 * (a_grid.grid[self.c.UCOMP,i] + a_grid.grid[self.c.UCOMP,i+1])/2

            P[i] = (a_grid.grid[self.c.UCOMP,i] * a_grid.grid[self.c.PCOMP, i+1]  + a_grid.grid[self.c.UCOMP,i+1] * a_grid.grid[self.c.PCOMP, i])/2

        for i in range(1, len(a_grid.x) - 1):
                       
            pressureEquilbriumFlux[self.c.RHOCOMP, i] = C1[i]

            pressureEquilbriumFlux[self.c.RHOCOMP + 1, i] = C2[i]

            pressureEquilbriumFlux[self.c.UCOMP, i] = Mu[i] + pi[i] 

            pressureEquilbriumFlux[self.c.PCOMP, i] =  K[i] + I[i] + P[i]

        for i in range(1, len(a_grid.x) - 1):

            leftWave = np.abs(a_grid.grid[self.c.UCOMP, i]) + np.sqrt(gammaBar[i] * a_grid.grid[self.c.PCOMP,i] / density[i])
            rightWave = np.abs(a_grid.grid[self.c.UCOMP, i+1]) + np.sqrt(gammaBar[i+1] * a_grid.grid[self.c.PCOMP,i+1] / density[i+1])
            sMaxR = max(leftWave, rightWave)


            leftWave = np.abs(a_grid.grid[self.c.UCOMP, i-1]) + np.sqrt(gammaBar[i-1] * a_grid.grid[self.c.PCOMP,i-1] / density[i-1])
            rightWave = np.abs(a_grid.grid[self.c.UCOMP, i]) + np.sqrt(gammaBar[i] * a_grid.grid[self.c.PCOMP,i] / density[i])
            sMaxL = max(leftWave, rightWave)
            for icomp in range(self.c.NUMQ):

                numFluxR[icomp, i] =  pressureEquilbriumFlux[icomp, i] - 0.5*sMaxR * (consU[icomp,i+1] - consU[icomp, i])
                numFluxL[icomp, i] = pressureEquilbriumFlux[icomp, i-1] - 0.5*sMaxL * (consU[icomp,i] - consU[icomp, i-1])

        return consU, numFluxR, numFluxL





    def keep_chandrashekar(self, a_grid, option=2):
        '''
        Kinetic Energy preserving and entropy stable finite volume schemes for compressible euler and navier stokes equations
        Chandrashekar

        Section 4.7 presents two options:
        1. Centered, kinetic energy preserving and entropy consistent numerical flux.
        2. Centered, kinetic energy preserving and entropy conservative numerical flux.
        '''

        
