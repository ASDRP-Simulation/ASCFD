from ascfd.euler import Euler
# from ascfd.reconstruct import weno5_reconstruction
from ascfd.constants import *
from ascfd.grid import Grid2D
from ascfd.constants import Constants

import numpy as np
import sys

class Flux:

    def __init__(self, a_constants: Constants, a_type: str):

        self.type = a_type
        self.c = a_constants

        self.euler = Euler(self.c)


        if self.type == "rusanov":
            self.flux_method = self.rusanov
        else:
            raise RuntimeError(f"Flux method not supported: {self.type}")


    def getFlux(self, a_grid):
        return self.flux_method(a_grid)


    def rusanov(self, a_grid):
        a_grid.assert_variable_type("prim")

        #get density 
        if self.c.system == "euler2D":
            density = a_grid.grid[self.c.RHOCOMP]
        else:
            raise RuntimeError("Density method needs to be implemented.")
            
        a = np.sqrt(self.c.gamma * a_grid.grid[self.c.PCOMP] / density)

        max_speed = np.max(np.abs(a_grid.grid[self.c.UCOMP]) + np.abs(a_grid.grid[self.c.VCOMP]) + a)

        #don't reconstruct cell faces just use FD methods
        U = a_grid.grid

        # conservative variables 
        consU = self.euler.prim_to_cons(U)

        # analytical flux
        fx, fy = self.euler.flux(U)

        #pointwise LF Flux
        numFluxX_plus = np.zeros_like(a_grid.grid)
        numFluxX_minus = np.zeros_like(a_grid.grid)
        numFluxY_plus = np.zeros_like(a_grid.grid)
        numFluxY_minus = np.zeros_like(a_grid.grid)

        # Computes LF Flux
        # Loop through all the cells except for the outermost ghost cells.
        for i in range(1, a_grid.Nx + 2*a_grid.Nghost - 1):
            for j in range(1, a_grid.Ny + 2*a_grid.Nghost - 1):
                # X-direction fluxes
                sMaxX_plus = max(
                    np.abs(a_grid.grid[self.c.UCOMP, i, j]) + a[i, j],
                    np.abs(a_grid.grid[self.c.UCOMP, i+1, j]) + a[i+1, j]
                )
                sMaxX_minus = max(
                    np.abs(a_grid.grid[self.c.UCOMP, i-1, j]) + a[i-1, j],
                    np.abs(a_grid.grid[self.c.UCOMP, i, j]) + a[i, j]
                )

                # Y-direction fluxes
                sMaxY_plus = max(
                    np.abs(a_grid.grid[self.c.VCOMP, i, j]) + a[i, j],
                    np.abs(a_grid.grid[self.c.VCOMP, i, j+1]) + a[i, j+1]
                )
                sMaxY_minus = max(
                    np.abs(a_grid.grid[self.c.VCOMP, i, j-1]) + a[i, j-1],
                    np.abs(a_grid.grid[self.c.VCOMP, i, j]) + a[i, j]
                )

                for icomp in range(self.c.NUMQ):
                    # Compute the Lax-Friedrichs flux for X-direction
                    numFluxX_plus[icomp, i, j] = 0.5*(fx[icomp, i+1, j] + fx[icomp, i, j]) - 0.5*sMaxX_plus * (consU[icomp, i+1, j] - consU[icomp, i, j])
                    numFluxX_minus[icomp, i, j] = 0.5*(fx[icomp, i, j] + fx[icomp, i-1, j]) - 0.5*sMaxX_minus * (consU[icomp, i, j] - consU[icomp, i-1, j])

                    # Compute the Lax-Friedrichs flux for Y-direction
                    numFluxY_plus[icomp, i, j] = 0.5*(fy[icomp, i, j+1] + fy[icomp, i, j]) - 0.5*sMaxY_plus * (consU[icomp, i, j+1] - consU[icomp, i, j])
                    numFluxY_minus[icomp, i, j] = 0.5*(fy[icomp, i, j] + fy[icomp, i, j-1]) - 0.5*sMaxY_minus * (consU[icomp, i, j] - consU[icomp, i, j-1])

        return consU, numFluxX_plus, numFluxX_minus, numFluxY_plus, numFluxY_minus
