from onedim.inputs import Inputs

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

        if inputs.system == "euler1D":
            #prim
            self.RHOCOMP = 0
            self.UCOMP = 1
            self.PCOMP = 2

            #cons
            self.RHOCOMP = 0
            self.MUCOMP = 1
            self.ECOMP = 2

            self.NUMQ = 3 
            self.variable_names = ["Density", "Velocity", "Pressure"]
            self.system = "euler1D"


            self.gammas = inputs.gammas #gammas
            self.mW = inputs.mW #molecular weight
            self.NS = 1 



        elif inputs.system == "euler2D":

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


        elif inputs.system == "euler1DNS2":
            
            self.NS = 2 #eventually don't hard code 2 species.

            #prim
            self.RHOCOMP = 0 #rhoY comps back to back
            self.UCOMP = self.NS
            self.PCOMP = self.NS + 1

            #cons
            self.RHOCOMP = 0 #rhoY 
            self.MUCOMP = self.NS
            self.ECOMP = self.NS + 1

            self.NUMQ = 4
            self.system = "euler1DNS2"
            self.variable_names = ["rhoY1", "rhoY2", "Velocity", "Pressure"]

            self.gammas = inputs.gammas #gammas
            self.mW = inputs.mW #molecular weight


        elif inputs.system == "4HRM":
            
            self.NS = 2 #eventually don't hard code 2 species.

            #prim
            self.RHOCOMP = 0 #rhoZ comps back to back
            self.UCOMP = self.NS
            self.PCOMP = self.NS + 1

            #cons
            self.RHOCOMP = 0 #rhoZ
            self.MUCOMP = self.NS
            self.ECOMP = self.NS + 1

            self.NUMQ = 4
            self.system = "euler1DNS2"
            self.variable_names = ["rhoZ1", "rhoZ2", "Velocity", "Pressure"]

            self.gammas = inputs.gammas #gammas
            self.mW = inputs.mW #molecular weight

        elif inputs.system == "euler2DNS2":
            raise RuntimeError(f"system in inputs file is not supported: {inputs.system}")            

        else:
            raise RuntimeError(f"system in inputs file is not supported: {inputs.system}")