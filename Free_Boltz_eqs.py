import numpy as np
import sys
from mpmath import *
from scipy.integrate import solve_ivp
import time
from Basic_definitions import *
from Quintuplet_definitions import *
from BoundStates_definitions import *
mp.dps = 25; mp.pretty = True


############################################################################

# Define Boltzmann equation without BS (Tree level)

############################################################################ 

def Boltzmann_Tree(DM, z, Y):

    prefactor = - DM.entropy(z)/( z * DM.Hubble(z) )

    abundances = Y**2 - DM.Yeq(z)**2

    return prefactor * DM.sv_production_NoSomm(z) * abundances


def do_scan_Boltzmann_Tree(M_DM_scan, print_results):

    Omega_Boltzmann_Tree = []

    for M_DM in M_DM_scan:
        
        DM = Quintuplet_DM(M_DM)

        bb = lambda z, Y: Boltzmann_Tree(DM, z, Y)

        solution = NDE_solver(bb, 0)

        omega_solution = Omega_DM(1.0, solution.y[0][-1], M_DM)

        Omega_Boltzmann_Tree.append([M_DM, omega_solution])

        if print_results == "y":
            print( [ M_DM, omega_solution ] )
        else: continue

    np.savetxt('./Results/Omega_Boltzmann_Tree.txt', Omega_Boltzmann_Tree)



############################################################################

# Define Boltzmann equation without BS (Sommerfeld with Hulthen potential)

############################################################################ 

def Boltzmann_Hulthen_free(DM, z, Y):

    prefactor = - DM.entropy(z)/( z * DM.Hubble(z) )

    abundances = Y**2 - DM.Yeq(z)**2

    return prefactor * DM.sv_production_Somm_Hulthen(z) * abundances


def do_scan_Boltzmann_Hulthen_free(M_DM_scan, print_results):

    Omega_Boltzmann_Hulthen_free = []

    for M_DM in M_DM_scan:
        
        DM = Quintuplet_DM(M_DM)

        bb = lambda z, Y: Boltzmann_Hulthen_free(DM, z, Y)

        solution = NDE_solver(bb, 0)

        omega_solution = Omega_DM(1.0, solution.y[0][-1], M_DM)

        Omega_Boltzmann_Hulthen_free.append([M_DM, omega_solution])

        if print_results == "y":
            print( [ M_DM, omega_solution ] )
        else: continue

    np.savetxt('./Results/Omega_Boltzmann_Hulthen_free.txt', Omega_Boltzmann_Hulthen_free)