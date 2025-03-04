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

# Define effective terms in Boltzmann equation

############################################################################ 


def BS_eff(BS, M, z):

    first_term = 1/BS.bsf(z)
    second_term_1 = gx**2 * sigma0_prime(M) * M**2/( 2 * BS.gI * BS.gamma_ann() )
    second_term_2 = ( 1 / (4*pi*z) )**(3/2)
    second_term_3 = np.exp( - z * BS.binding_energy_BS(z) )

    return sigma0_prime(M) * 1/(first_term + second_term_1 * second_term_2 * second_term_3)


def Boltzmann_Hulthen_effective(DM, BS_list, z, Y):

    M = DM.M

    prefactor = - DM.entropy(z)/( z * DM.Hubble(z) )

    abundances = Y**2 - DM.Yeq(z)**2

    sv_eff = DM.sv_production_Somm_Hulthen(z) + np.sum([BS_eff(BS, M, z) for BS in BS_list])

    return prefactor * sv_eff * abundances


############################################################################

# Define Boltzmann equation with 1s BS (effective equation)

############################################################################ 

def do_scan_Boltzmann_Hulthen_1s_effective(M_DM_scan, print_results):

    Omega_Boltzmann_Hulthen_1s_effective = []

    for M_DM in M_DM_scan:
        
        DM = Quintuplet_DM(M_DM)

        # # Variable order is: (Mass, gI, nE, l, I, nS, BS function, group factor constant)

        BS_1s1 = Quintuplet_BS(M_DM, gI_list[0], nE_list[0], l_list[0], Isp_list[0], nS_Quint, BS_func_list[0], gf_list[0])
        BS_1s3 = Quintuplet_BS(M_DM, gI_list[1], nE_list[1], l_list[1], Isp_list[1], nS_Quint, BS_func_list[1], gf_list[1])
        BS_1s5 = Quintuplet_BS(M_DM, gI_list[2], nE_list[2], l_list[2], Isp_list[2], nS_Quint, BS_func_list[2], gf_list[2])

        BS_1s_list = [BS_1s1, BS_1s3, BS_1s5]

        bb = lambda z, Y: Boltzmann_Hulthen_effective(DM, BS_1s_list, z, Y)

        solution = NDE_solver(bb, 0)
        
        omega_solution = Omega_DM(1.0, solution.y[0][-1], M_DM)

        Omega_Boltzmann_Hulthen_1s_effective.append([M_DM, omega_solution])

        if print_results == "y":
            print( [ M_DM, omega_solution ] )
        else: continue

    np.savetxt('./Results/Omega_Boltzmann_Hulthen_1s_effective.txt', Omega_Boltzmann_Hulthen_1s_effective)



############################################################################

# Define Boltzmann equation with 1s, 2s, 2p BS (effective equation)

############################################################################ 



def do_scan_Boltzmann_Hulthen_1s2s2p_effective(M_DM_scan, print_results):

    Omega_Boltzmann_Hulthen_1s2s2p_effective = []

    for M_DM in M_DM_scan:
        
        DM = Quintuplet_DM(M_DM)

        # # Variable order is: (Mass, gI, nE, l, I, nS, BS function, group factor constant)

        BS_1s1 = Quintuplet_BS(M_DM, gI_list[0], nE_list[0], l_list[0], Isp_list[0], nS_Quint, BS_func_list[0], gf_list[0])
        BS_1s3 = Quintuplet_BS(M_DM, gI_list[1], nE_list[1], l_list[1], Isp_list[1], nS_Quint, BS_func_list[1], gf_list[1])
        BS_1s5 = Quintuplet_BS(M_DM, gI_list[2], nE_list[2], l_list[2], Isp_list[2], nS_Quint, BS_func_list[2], gf_list[2])

        BS_2s1 = Quintuplet_BS(M_DM, gI_list[3], nE_list[3], l_list[3], Isp_list[3], nS_Quint, BS_func_list[3], gf_list[3])
        BS_2s3 = Quintuplet_BS(M_DM, gI_list[4], nE_list[4], l_list[4], Isp_list[4], nS_Quint, BS_func_list[4], gf_list[4])
        BS_2s5 = Quintuplet_BS(M_DM, gI_list[5], nE_list[5], l_list[5], Isp_list[5], nS_Quint, BS_func_list[5], gf_list[5])

        BS_2p1 = Quintuplet_BS(M_DM, gI_list[6], nE_list[6], l_list[6], Isp_list[6], nS_Quint, BS_func_list[6], gf_list[6])
        BS_2p3 = Quintuplet_BS(M_DM, gI_list[7], nE_list[7], l_list[7], Isp_list[7], nS_Quint, BS_func_list[7], gf_list[7])
        BS_2p5 = Quintuplet_BS(M_DM, gI_list[8], nE_list[8], l_list[8], Isp_list[8], nS_Quint, BS_func_list[8], gf_list[8])

        BS_list = [BS_1s1, BS_1s3, BS_1s5, BS_2s1, BS_2s3, BS_2s5, BS_2p1, BS_2p3, BS_2p5]

        bb = lambda z, Y: Boltzmann_Hulthen_effective(DM, BS_list, z, Y)

        solution = NDE_solver(bb, 0)
        
        omega_solution = Omega_DM(1.0, solution.y[0][-1], M_DM)

        Omega_Boltzmann_Hulthen_1s2s2p_effective.append([M_DM, omega_solution])

        if print_results == "y":
            print( [ M_DM, omega_solution ] )
        else: continue

    np.savetxt('./Results/Omega_Boltzmann_Hulthen_1s2s2p_effective.txt', Omega_Boltzmann_Hulthen_1s2s2p_effective)