import numpy as np
import matplotlib.pyplot as plt
from ascfd.constants import Constants


class Grid1D:
    def __init__(self, xlim, Nx, Nghost, num_vars):

        self.num_vars = num_vars
        
        self.xlim = xlim
        self.Nx = Nx
        self.Nghost = Nghost
        self.num_vars = num_vars  # Number of variables (e.g., 3 for Euler equations)
        self.dx = (xlim[1] - xlim[0]) / (Nx - 1)
        self.x = np.linspace( 
            xlim[0] - self.dx * Nghost, xlim[1] + self.dx * Nghost, Nx + 2 * Nghost
        )
        self.ndim = 1
        # Internal x (no ghost cells)
        self.x_int = self.x[self.Nghost : self.Nx + self.Nghost]
        # Initialize the grid for multiple variables
        self.grid = np.zeros((num_vars, Nx + 2 * Nghost))
        self.variables = "prim"

        self.shock_mask = np.zeros_like(self.x)

    def fill_grid(self, f):
        for var in range(self.num_vars):
            self.grid[var, self.Nghost : self.Nx + self.Nghost] = f(self.x_int, var)

    def plot(self):
        fig, axs = plt.subplots(self.num_vars, 1, figsize=(10, self.num_vars * 3))
        if self.num_vars == 1:
            axs = [axs]  # Make it iterable
        for var in range(self.num_vars):
            axs[var].plot(self.x, self.grid[var, :], ".")
            axs[var].set_title(self.c.variable_names[var])
        plt.tight_layout()
        plt.show()


    def return_internal_grid(self):
        return self.grid[:, self.Nghost : self.Nx + self.Nghost]

    def return_grid(self):
        return self.grid

    def scratch_internal(self):
        return np.zeros_like(self.grid[:, self.Nghost : self.Nx + self.Nghost])

    def scratch(self):
        return np.zeros_like(self.grid)

    def transform(self, f, name):
        self.grid = f(self.grid)
        self.variables = name

    def set_internal(self, U_new):
        self.grid[:, self.Nghost : self.Nx + self.Nghost] = U_new

    def set(self, U_new):
        assert np.shape(U_new) == np.shape(self.grid)
        self.grid = U_new

    def set_grid(self, U_new):
        self.grid = U_new

    def internal_bounds(self):
        return (self.Nghost, self.Nx + self.Nghost)

    def assert_variable_type(self, str):
        assert str == self.variables

    def check_grid(self, constants, prim=False, cons=False):
        for i in range(len(self.grid[0])):

            if prim:
                if self.grid[constants.PCOMP, i] <= 0:
                    print("Negative Pressure - Bad cell: ", i)
                    assert False

                if self.grid[constants.RHOCOMP, i] <= 0:
                    print("Negative Pressure - Bad cell: ", i)
                    assert False

            if cons:
                if self.grid[constants.ECOMP, i] <= 0:
                    print("Negative Energy - Bad cell: ", i)
                    assert False


            for icomp in range(constants.NUMQ):
                if np.isnan(self.grid[icomp, i]):
                    print("Nan cell : ", i)
                    assert False



class Grid2D:
    def __init__(self, xlim, ylim, Nx, Ny, Nghost, num_vars):
        self.xlim = xlim
        self.ylim = ylim
        self.Nx = Nx
        self.Ny = Ny
        self.Nghost = Nghost
        self.num_vars = num_vars  
        self.dx = (xlim[1] - xlim[0]) / (Nx - 1)
        self.dy = (ylim[1] - ylim[0]) / (Ny - 1)
        self.x = np.linspace(
            xlim[0] - self.dx * Nghost, xlim[1] + self.dx * Nghost, Nx + 2 * Nghost
        )
        self.y = np.linspace(
            ylim[0] - self.dy * Nghost, ylim[1] + self.dy * Nghost, Ny + 2 * Nghost
        )

        self.meshX, self.meshY = np.meshgrid(self.x, self.y)

        self.ndim = 2
        # Internal x (no ghost cells)
        self.x_int = self.x[self.Nghost : self.Nx + self.Nghost]
        self.y_int = self.y[self.Nghost : self.Ny + self.Nghost]

        self.meshIntX, self.meshIntY = np.meshgrid(self.x_int, self.y_int)

        # Initialize the grid for multiple variables
        self.grid = np.zeros((num_vars, Nx + 2 * Nghost, Ny + 2 * Nghost))
        self.variables = "prim"

    def fill_grid(self, f):
        for var in range(self.num_vars):
            # we can jsut fill ghost cells too, i don't think it matters we're just going to override anyways
            self.grid[var] = f(self.meshX, self.meshY, var)

    def plot(self):
        # Create a figure with 4 subplots arranged in 2x2
        fig, axs = plt.subplots(2, 2, figsize=(12, 10))
        axs = axs.flatten()  # Flatten the 2x2 grid into a 1D array of axes

        # Titles for each subplot corresponding to each variable
        titles = ["Density", "UVel", "VVel", "Pressure"]

        # Loop over the number of variables and plot each one using the 'plasma' colormap
        for i in range(self.num_vars):
            pcm = axs[i].pcolormesh(
                self.meshX, self.meshY, self.grid[i].T, cmap="inferno", shading="auto"
            )
            fig.colorbar(
                pcm, ax=axs[i], orientation="vertical"
            )  # Add a color bar to each subplot
            axs[i].set_title(titles[i])
            axs[i].set_xlabel("X")
            axs[i].set_ylabel("Y")

        plt.tight_layout()  # Adjust layout to prevent overlap
        plt.show()  # Display the plot

    def apply_periodic_bcs(self):
        for var in range(self.num_vars):
            for i in range(self.Nghost):
                self.grid[var, i, :] = self.grid[
                    var, self.Nx + self.Nghost - (self.Nghost - i), :
                ]
                self.grid[var, self.Nx + self.Nghost + i, :] = self.grid[
                    var, self.Nghost + i, :
                ]

            for j in range(self.Nghost):
                self.grid[var, :, j] = self.grid[
                    var, :, self.Ny + self.Nghost - (self.Nghost - j)
                ]
                self.grid[var, :, self.Ny + self.Nghost + j] = self.grid[
                    var, :, self.Nghost + j
                ]

    def apply_zero_gradient_bcs(self):
        for var in range(self.num_vars):
            for ighost in range(self.Nghost):
                self.grid[var, ighost, :] = self.grid[var, self.Nghost, :]
                self.grid[var, self.Nx + self.Nghost + ighost, :] = self.grid[
                    var, self.Nx + self.Nghost - 1, :
                ]

            for jghost in range(self.Nghost):
                self.grid[var, :, jghost] = self.grid[var, :, self.Nghost]
                self.grid[var, :, self.Ny + self.Nghost + jghost] = self.grid[
                    var, :, self.Ny + self.Nghost - 1
                ]

    def return_internal_grid(self):
        pass

    def return_grid(self):
        pass

    def scratch_internal(self):
        pass

    def scratch(self):
        pass

    def transform(self, f, name):
        pass

    def set_internal(self, U_new):
        pass

    def set(self, U_new):
        pass

    def set_grid(self, U_new):
        pass

    def internal_bounds(self):
        pass

    def assert_variable_type(self, str):
        pass
