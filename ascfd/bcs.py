
from ascfd.grid import Grid1D
from typing import Callable

class BoundaryConditions:
    
    def __init__(self, grid: Grid1D, type_lo: str, type_hi: str) -> None:

        self.type_lo = type_lo
        self.type_hi = type_hi

        #mutable objects, such as Grid1D are passed by reference in python
        self.grid = grid

        if self.type_lo == "dirichlet":
            self.f_lo = self.dirichlet_lo
        elif self.type_lo == "neumann":
            self.f_lo = self.neumann_lo
        elif self.type_lo == "periodic":
            self.f_lo = self.periodic_lo
        else:
            raise RuntimeError(f"BC Lo Type not supported: {self.type_lo}")

        if self.type_hi == "dirichlet":
            self.f_hi = self.dirichlet_hi
        elif self.type_hi == "neumann":
            self.f_hi = self.neumann_hi
        elif self.type_hi == "periodic":
            self.f_hi = self.periodic_hi
        else:
            raise RuntimeError(f"BC Hi Type not supported: {self.type_hi}")


    def apply_bcs(self) -> None:
        self.apply_lo()
        self.apply_hi()

    def apply_lo(self) -> None:    
        self.f_lo(self.grid)
        
    def apply_hi(self) -> None:
        self.f_hi(self.grid)


    @staticmethod
    def dirichlet_lo(grid: Grid1D) -> None:
        for var in range(grid.num_vars):
            grid.grid[var, :grid.Nghost] = 0  

    @staticmethod
    def dirichlet_hi(grid: Grid1D) -> None:
        for var in range(grid.num_vars):
            grid.grid[var, -grid.Nghost:] = 0 

    @staticmethod
    def neumann_lo(grid: Grid1D) -> None:
        for var in range(grid.num_vars):
            grid.grid[var, :grid.Nghost] = grid.grid[var, grid.Nghost] 

    @staticmethod
    def neumann_hi(grid: Grid1D) -> None:
        for var in range(grid.num_vars):
            grid.grid[var, -grid.Nghost:] = grid.grid[var, -grid.Nghost-1] 

    @staticmethod
    def periodic_lo(grid: Grid1D) -> None:
        for var in range(grid.num_vars):
            grid.grid[var, :grid.Nghost] = grid.grid[var, -2*grid.Nghost:-grid.Nghost]

    @staticmethod
    def periodic_hi(grid: Grid1D) -> None:
        for var in range(grid.num_vars):
            grid.grid[var, -grid.Nghost:] = grid.grid[var, grid.Nghost:2*grid.Nghost]  
