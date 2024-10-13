from ascfd.inputs import Inputs

# primative
RHOCOMP = 0
UCOMP = 1
PCOMP = 2

# conservative
RHOCOMP = 0
MUCOMP = 1
ECOMP = 2


NUMQ = 3

variable_names = ["Density", "Velocity", "Pressure"]


# Fluid
gamma = 1.4



class Constants:

    def __init__(self, inputs: Inputs):

        self.gamma = 1.4
        self.system = inputs.system

        if inputs.system == "euler2D":

            #prim
            self.RHOCOMP = 0
            self.UCOMP = 1
            self.VCOMP = 2
            self.PCOMP = 3
 
            #cons
            self.RHOCOMP = 0
            self.MUCOMP = 1
            self.MVCOMP = 2
            self.ECOMP = 3

            self.NUMQ = 4

            self.system = "euler2D"

            self.gammas = inputs.gammas #gammas
            self.mW = inputs.mW #molecular weight

            self.NS = 1 

        else:
            raise RuntimeError(f"system in inputs file is not supported: {inputs.system}")