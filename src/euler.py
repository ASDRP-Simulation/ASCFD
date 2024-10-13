import numpy as np
from onedim.constants import Constants

class Euler:

    def __init__(self, a_constants : Constants):
    
        self.c = a_constants


    # def get_gamma_bar(self, a_Y):
        
    #     R_u = 8.314# kJ / (kmol. K)
        
    #     cps = np.zeros(self.c.NS)
    #     cvs = np.zeros(self.c.NS)

    #     for i in range(self.c.NS):
    #         #mixing rule eq 4 in conservative presesure equilbirum
    #         cps[i] = (self.c.gammas[i] / (self.c.gammas[i] - 1)) * (R_u/ self.c.mW[i])
    #         cvs[i] = (1 / (self.c.gammas[i] - 1)) * (R_u/ self.c.mW[i])
        
    #     top = 0
    #     bottom = 0
    #     for i in range(self.c.NS):
    #         top += cps[i] * a_Y[i]
    #         bottom += cvs[i] * a_Y[i]


    #     return top/bottom
        




    def prim_to_cons(self, a_prim):
        cons = np.zeros_like(a_prim)

        if self.c.system == "euler1D":

            cons[self.c.RHOCOMP] = a_prim[self.c.RHOCOMP]
            cons[self.c.MUCOMP] = a_prim[self.c.RHOCOMP] * a_prim[self.c.UCOMP]
            E = a_prim[self.c.PCOMP] / ((self.c.gamma - 1) * a_prim[self.c.RHOCOMP]) + 0.5 * a_prim[self.c.UCOMP] ** 2
            cons[self.c.ECOMP] = E * a_prim[self.c.RHOCOMP]
        
        elif self.c.system == "euler1DNS2":

            totalDensity = np.zeros_like(cons[0])
            for k in range(self.c.NS):
                cons[k] = a_prim[k]                
                totalDensity += cons[k]
            
            cons[self.c.MUCOMP] = totalDensity*a_prim[self.c.UCOMP]
            # E = a_prim[self.c.PCOMP] / ((self.c.gamma - 1) * totalDensity) + 0.5 * a_prim[self.c.UCOMP] ** 2
            # cons[self.c.ECOMP] = E * totalDensity

            # rhoY / rho = Y
            Ys = [ a_prim[self.c.RHOCOMP] / totalDensity, a_prim[self.c.RHOCOMP+1] / totalDensity ]

            # gammaBar = self.get_gamma_bar(Ys)


            #eq 29
            mBar = totalDensity / ( a_prim[self.c.RHOCOMP] / self.c.mW[0] + a_prim[self.c.RHOCOMP + 1] / self.c.mW[1])

            #eq 4
            sum = 0
            for i in range(self.c.NS):

                #Y = rhoY/rho
                Y_i = a_prim[self.c.RHOCOMP + i]/totalDensity

                sum += (1 / (1 - self.c.gammas[i])) * (Y_i / self.c.mW[i])

            G = mBar * sum
            gammaBar = 1/G + 1



                
            e = a_prim[self.c.PCOMP] / ((gammaBar - 1) * totalDensity)
            E = totalDensity * e + 0.5 * totalDensity * a_prim[self.c.UCOMP] ** 2
            cons[self.c.ECOMP] = E


        elif self.c.system == "4HRM":

            totalDensity = np.zeros_like(cons[0])
            for k in range(self.c.NS):
                cons[k] = a_prim[k]                
                totalDensity += cons[k]
            
            cons[self.c.MUCOMP] = totalDensity*a_prim[self.c.UCOMP]
            # E = a_prim[self.c.PCOMP] / ((self.c.gamma - 1) * totalDensity) + 0.5 * a_prim[self.c.UCOMP] ** 2
            # cons[self.c.ECOMP] = E * totalDensity

            ## rhoZ / rho = Z
            #Ys = [ a_prim[self.c.RHOCOMP] / totalDensity, a_prim[self.c.RHOCOMP+1] / totalDensity ]

            #rho_k alpha_k / rho = Y_k
            Y_1 = [a_prim[self.c.RHOCOMP] / totalDensity, a_prim[self.c.RHOCOMP] / totalDensity]
            Y_2 = [a_prim[self.c.RHOCOMP] / totalDensity, a_prim[self.c.RHOCOMP] / totalDensity]

            Ys = [Y_1, Y_2 ]


            

            # gammaBar = self.get_gamma_bar(Ys)


            #eq 29
            mBar = totalDensity / ( a_prim[self.c.RHOCOMP] / self.c.mW[0] + a_prim[self.c.RHOCOMP + 1] / self.c.mW[1])

            #eq 4
            sum = 0
            for i in range(self.c.NS):

                #Y = rhoY/rho
                Y_i = a_prim[self.c.RHOCOMP + i]/totalDensity

                sum += (1 / (1 - self.c.gammas[i])) * (Y_i / self.c.mW[i])

            G = mBar * sum
            gammaBar = 1/G + 1



                
            e = a_prim[self.c.PCOMP] / ((gammaBar - 1) * totalDensity)
            E = totalDensity * e + 0.5 * totalDensity * a_prim[self.c.UCOMP] ** 2
            cons[self.c.ECOMP] = E

            

        else:
            raise RuntimeError(f"System not supported: {self.c.system}")    
        

        return cons




    def cons_to_prim(self, a_cons):
        prim = np.zeros_like(a_cons)


        if self.c.system == "euler1D":

            prim[self.c.RHOCOMP] = a_cons[self.c.RHOCOMP]
            prim[self.c.UCOMP] = a_cons[self.c.MUCOMP] / a_cons[self.c.RHOCOMP]
            prim[self.c.PCOMP] = (self.c.gamma - 1) * (
                a_cons[self.c.ECOMP] - 0.5 * a_cons[self.c.RHOCOMP] * prim[self.c.UCOMP] ** 2
            )
        elif self.c.system == "euler1DNS2":

            totalDensity = np.zeros_like(a_cons[0])
            for k in range(self.c.NS):
                prim[k] = a_cons[k]                
                totalDensity += a_cons[k]
            prim[self.c.UCOMP] = a_cons[self.c.MUCOMP] / totalDensity

            # rhoY / rho = Y
            Ys = [ a_cons[self.c.RHOCOMP] / totalDensity, a_cons[self.c.RHOCOMP+1] / totalDensity ]

            # gammaBar = self.get_gamma_bar(Ys)

            #eq 29
            mBar = totalDensity / ( a_cons[self.c.RHOCOMP] / self.c.mW[0] + a_cons[self.c.RHOCOMP + 1] / self.c.mW[1])

            #eq 4
            sum = 0
            for i in range(self.c.NS):

                #Y = rhoY/rho
                Y_i = a_cons[self.c.RHOCOMP + i]/totalDensity

                sum += (1 / (1 - self.c.gammas[i])) * (Y_i / self.c.mW[i])

            G = mBar * sum
            gammaBar = 1/G + 1


            prim[self.c.PCOMP] = (gammaBar - 1) * \
                                 (a_cons[self.c.ECOMP] - 0.5 * totalDensity * prim[self.c.UCOMP] ** 2)


            
        else:
            raise RuntimeError(f"System not supported: {self.c.system}")



        return prim



    def flux(self, a_prim):
        """
        Compute the flux for the Euler equations
        """

        flux = np.zeros_like(a_prim)

        if self.c.system == "euler1D":

            flux[self.c.RHOCOMP, :] = a_prim[self.c.RHOCOMP] * a_prim[self.c.UCOMP]
            flux[self.c.UCOMP, :] = a_prim[self.c.RHOCOMP] * a_prim[self.c.UCOMP] ** 2 + a_prim[self.c.PCOMP]

            e = a_prim[self.c.PCOMP] / ((self.c.gamma - 1) * a_prim[self.c.RHOCOMP])
            E = a_prim[self.c.RHOCOMP] * e + 0.5 * a_prim[self.c.RHOCOMP] * a_prim[self.c.UCOMP] ** 2
            flux[self.c.PCOMP, :] = (E + a_prim[self.c.PCOMP]) * a_prim[self.c.UCOMP]

        elif self.c.system == "euler1DNS2":

            totalDensity = np.zeros_like(a_prim[0])
            for k in range(self.c.NS):
                flux[k] = a_prim[k] * a_prim[self.c.UCOMP]          
                totalDensity += a_prim[k]
            
            
            flux[self.c.UCOMP] = totalDensity* (a_prim[self.c.UCOMP])**2 + a_prim[self.c.PCOMP]

            # rhoY / rho = Y
            Ys = [ a_prim[self.c.RHOCOMP] / totalDensity, a_prim[self.c.RHOCOMP+1] / totalDensity ]

            #eq 29
            mBar = totalDensity / ( a_prim[self.c.RHOCOMP] / self.c.mW[0] + a_prim[self.c.RHOCOMP + 1] / self.c.mW[1])

            #eq 4
            sum = 0
            for i in range(self.c.NS):

                #Y = rhoY/rho
                Y_i = a_prim[self.c.RHOCOMP + i]/totalDensity

                sum += (1 / (1 - self.c.gammas[i])) * (Y_i / self.c.mW[i])

            G = mBar * sum
            gammaBar = 1/G + 1



            
            e = a_prim[self.c.PCOMP] / ((gammaBar - 1) * totalDensity)
            E = totalDensity * e + 0.5 * totalDensity * a_prim[self.c.UCOMP] ** 2

            flux[self.c.PCOMP, :] = (E + a_prim[self.c.PCOMP]) * a_prim[self.c.UCOMP]



        return flux


 
    def get_max_speed(self, a_grid):
        if a_grid.variables == "prim":
            return np.max(a_grid.grid[self.c.UCOMP])
        elif a_grid.variables == "cons":
            return np.max(a_grid.grid[self.c.MUCOMP] / a_grid.grid[self.c.RHOCOMP])
        else:
            print("unsupported")
            exit()


    # def prim_to_char(self, prim_vars):
    #     num_points = prim_vars.shape[1]

    #     # Initialize the array for characteristic variables
    #     char_vars = np.zeros_like(prim_vars)

    #     # Initialize a list or 3D array to store the right eigenvectors for each point
    #     right_eigenvectors_list = []

    #     for i in range(num_points):
    #         rho = prim_vars[self.c.RHOCOMP, i]
    #         u = prim_vars[self.c.UCOMP, i]
    #         p = prim_vars[self.c.PCOMP, i]

    #         A = np.array(
    #             [
    #                 [0, 1, 0],
    #                 [-(u**2) + self.c.gamma * p / rho, 2 * u, self.c.gamma - 1],
    #                 [
    #                     -(self.c.gamma - 1) * u**3 + self.c.gamma * u * p / rho,
    #                     self.c.gamma * u**2 - 1.5 * (self.c.gamma - 1) * u**2,
    #                     self.c.gamma * u,
    #                 ],
    #             ]
    #         )

    #         eigenvalues, right_eigenvectors = np.linalg.eig(A)

    #         # Store the right eigenvectors for each point
    #         right_eigenvectors_list.append(right_eigenvectors)

    #         # Inverse of the right eigenvectors matrix
    #         R_inv = np.linalg.inv(right_eigenvectors)

    #         U = np.array([rho, u, p])
    #         char_vars[:, i] = R_inv @ U

    #     return char_vars, np.array(right_eigenvectors_list)


    # def char_to_prim(self, char_vars, right_eigenvectors_list):
    #     # Assuming char_vars is a 2D array with shape (3, N) where N is the number of spatial points
    #     num_points = char_vars.shape[1]

    #     # Initialize the array for primitive variables
    #     prim_vars = np.zeros_like(char_vars)

    #     for i in range(num_points):
    #         # Get the right eigenvectors for the current point
    #         right_eigenvectors = right_eigenvectors_list[i]

    #         # Transform characteristic variables back to primitive variables at each point
    #         prim_vars[:, i] = right_eigenvectors @ char_vars[:, i]

    #     return prim_vars
