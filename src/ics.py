from onedim.constants import *
import numpy as np


def sod_shock_tube(a_x, a_var):
    ics = np.zeros_like(a_x)

    for i in range(len(a_x)):
        if a_var == RHOCOMP:
            if a_x[i] < 0.5:
                ics[i] = 1
            else:
                ics[i] = 0.125

        elif a_var == PCOMP:
            if a_x[i] < 0.5:
                ics[i] = 1
            else:
                ics[i] = 0.1

        elif a_var == UCOMP:
            if a_x[i] < 0.5:
                ics[i] = 0
            else:
                ics[i] = 0

        else:
            print("Unexpected Variable")
            exit()
    return ics


def shu_osher_shock_tube(a_x, a_var):

    #https://arxiv.org/pdf/1711.11288
    ics = np.zeros_like(a_x)

    shock_pos = -4

    for i in range(len(a_x)):
        if a_var == RHOCOMP:
            if a_x[i] < shock_pos:
                ics[i] = 27/7
            else:
                ics[i] = 1 + 0.2 * np.sin(5*np.pi*a_x[i])

        elif a_var == UCOMP:
            if a_x[i] < shock_pos:
                ics[i] = 4*np.sqrt(35)/9
            else:
                ics[i] = 0.0


        elif a_var == PCOMP:
            if a_x[i] < shock_pos:
                ics[i] = 31/3
            else:
                ics[i] = 1.0

        else:
            print("Unexpected Variable")
            exit()
    return ics


def lax_problem(a_x, a_var):

    #https://arxiv.org/pdf/1711.11288
    ics = np.zeros_like(a_x)

    shock_pos = 0

    for i in range(len(a_x)):
        if a_var == RHOCOMP:
            if a_x[i] < shock_pos:
                ics[i] = 0.445
            else:
                ics[i] = 0.5

        elif a_var == UCOMP:
            if a_x[i] < shock_pos:
                ics[i] = 0.698
            else:
                ics[i] = 0.0


        elif a_var == PCOMP:
            if a_x[i] < shock_pos:
                ics[i] = 3.528
            else:
                ics[i] = 0.571

        else:
            print("Unexpected Variable")
            exit()
    return ics


def NS2_smooth_interfaces(a_x, a_var):


    # Fully conservative and pressure-equilibrium preserving
    #scheme for compressible multi-component flows
    #Yuji Fujiwara, Yoshiharu Tamaki, Soshi Kawai∗
    #4.1. 1D inviscid smooth interfaces advection problem in pressure and velocity equilibriums
    
    #parameters used in paper
    xc = 0.5
    rc = 0.25
    (w1, w2) = (0.6, 0.2)
    k = 20


    ics = np.zeros_like(a_x)

    for i in range(len(a_x)):
        if a_var == 0: #rhoY1
            ics[i] = (w1/2) * (1 - np.tanh(k*(np.abs(a_x[i] - xc)  - rc)))

        elif a_var == 1: #rhoY2
            ics[i] = (w2/2) * (1 + np.tanh(k*(np.abs(a_x[i] - xc)  - rc)))


        elif a_var == 2: #u
            ics[i] = 1
        
        elif a_var == 3: #p
            ics[i] = 0.9

        else:
            print("Unexpected Variable")
            exit()
    return ics



    return ics

