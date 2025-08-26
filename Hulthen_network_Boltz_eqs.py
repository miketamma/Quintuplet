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

# Define ratios for the network cases

############################################################################

# # def Yratio(BS, M, z):

# #     prefactor = ( pi**(5/2) ) * 32/45

# #     prefactor_g = BS.gI * gs_star(M/z)/( gx**2 )

# #     factor_z = ( np.exp( z * BS.binding_energy_BS(z) ) )/( z**(3/2) )

# #     return prefactor * prefactor_g * factor_z

# def Yratio(BS, M, z):

#     prefactor = ( pi**(5/2) ) * 32/45

#     prefactor_g = BS.gI * gs_star(M/z)/( gx**2 )

#     factor_z = ( np.exp( - z * BS.binding_energy_BS(z) ) ) * ( z**(3/2) )

#     return factor_z/( prefactor * prefactor_g )

def Yratio(BS, M, z):

    prefactor = ( pi**(7/2) ) * 16/45

    prefactor_g = BS.gI * gs_star(M/z)/( gx**2 )

    factor_z = ( np.exp( z * BS.binding_energy_BS(z) ) )/( z**(3/2) )

    return prefactor * prefactor_g * factor_z

def InverseYratio(BS, M, z):

    # prefactor = ( pi**(5/2) ) * 32/45

    # prefactor_g = BS.gI * gs_star(M/z)/( gx**2 )

    # factor_z = ( np.exp( - z * BS.binding_energy_BS(z) ) ) * ( z**(3/2) )

    return 1/Yratio(BS, M, z)

def Yratio_NoExp(BS, M, z):

    prefactor = ( pi**(7/2) ) * 16/45

    prefactor_g = BS.gI * gs_star(M/z)/( gx**2 )

    factor_z = 1/( z**(3/2) )

    return prefactor * prefactor_g * factor_z




############################################################################

# Define Boltzmann equation Network with 1 BS

############################################################################ 


def Boltzmann_Hulthen_Network_1s1(DM, BS, z, Y):


    Y_DM = Y[0]
    Y_1s1 = Y[1]

    ### Define DM mass ###

    M = DM.M

    ### Boltzmann equation for DM ###

    prefactor_DM = - DM.entropy(z)/( z * DM.Hubble(z) )
    prefactor_I = 1/( z * DM.Hubble(z) )

    dYdz_DM_1 = ( Y_DM**2 - DM.Yeq(z)**2 ) * DM.sv_production_Somm_Hulthen(z)
    dYdz_DM_2 = sigma0_prime(M) * BS.bsf(z) * ( Y_DM**2 - Y_1s1 * InverseYratio(BS, M, z) )

    dYdz_DM = prefactor_DM * ( dYdz_DM_1 + dYdz_DM_2 )

    ### Boltzmann equation for BS ###

    prefactor_I = 1/( z * DM.Hubble(z) )

    dYdz_1s1_1 = BS.Gamma_break_NoExp(z) * Yratio_NoExp(BS, M, z) * Y_DM**2
    dYdz_1s1_2 = - BS.Gamma_break(z) * Y_1s1
    dYdz_1s1_3 = BS.Gamma_ann_Hulthen_TA(z) * ( BS.Y_BS_eq(z) - Y_1s1 )

    dYdz_1s1 = prefactor_I * ( dYdz_1s1_1 + dYdz_1s1_2 + dYdz_1s1_3)

    return [dYdz_DM, dYdz_1s1]


def do_scan_Boltzmann_Hulthen_Network_1s1(M_DM_scan, print_results):

    Omega_Boltzmann_Hulthen_Network_1s1 = []

    for M_DM in M_DM_scan:
        
        DM = Quintuplet_DM(M_DM)

        # # Variable order is: (Mass, gI, nE, l, I, nS, BS function, group factor constant)

        BS_1s1 = Quintuplet_BS(M_DM, gI_list[0], nE_list[0], l_list[0], Isp_list[0], nS_Quint, BS_func_list[0], gf_list[0])

        bb = lambda z, Y: Boltzmann_Hulthen_Network_1s1(DM, BS_1s1, z, Y)

        solution = NDE_solver(bb, 1)
        
        omega_solution = Omega_DM(1.0, solution.y[0][-1], M_DM)

        Omega_Boltzmann_Hulthen_Network_1s1.append([M_DM, omega_solution])

        if print_results == "y":
            print( [ M_DM, omega_solution ] )
        else: continue

    np.savetxt('./Results/Omega_Boltzmann_Hulthen_Network_1s1.txt', Omega_Boltzmann_Hulthen_Network_1s1)


############################################################################

# Define Boltzmann equation Network with all 1s BS

############################################################################ 


def Boltzmann_Hulthen_Network_1s(DM, BS_list, z, Y):


    Y_DM = Y[0]
    Y_1s1 = Y[1]
    Y_1s3 = Y[2]
    Y_1s5 = Y[3]
    Y_1s = [Y_1s1, Y_1s3, Y_1s5]

    ### Define DM mass ###

    M = DM.M

    ### Boltzmann equation for DM ###

    prefactor_DM = - DM.entropy(z)/( z * DM.Hubble(z) )

    dYdz_DM_1 = ( Y_DM**2 - DM.Yeq(z)**2 ) * DM.sv_production_Somm_Hulthen(z)
    dYdz_DM_2 = sigma0_prime(M) * np.sum( [ BS_list[i].bsf(z) * ( Y_DM**2 - Y_1s[i] * InverseYratio(BS_list[i], M, z) ) for i in range(len(BS_list)) ] )

    dYdz_DM = prefactor_DM * ( dYdz_DM_1 + dYdz_DM_2 )

    ### Boltzmann equation for BS ###

    prefactor_I = 1/( z * DM.Hubble(z) )

    dYdz_1s = []

    for i in range(len(BS_list)):
        dYdz_1s_1 = BS_list[i].Gamma_break_NoExp(z) * Yratio_NoExp(BS_list[i], M, z) * Y_DM**2 
        dYdz_1s_2 = BS_list[i].Gamma_break(z) * Y_1s[i]  
        dYdz_1s_3 = BS_list[i].Gamma_ann_Hulthen_TA(z) * ( BS_list[i].Y_BS_eq(z) - Y_1s[i] )

        dYdz_1s.append( prefactor_I * ( dYdz_1s_1 - dYdz_1s_2 + dYdz_1s_3 ) )

    # dYdz_1s_1 = [ BS_list[i].Gamma_break_NoExp(z) * Yratio_NoExp(BS_list[i], M, z) * Y_DM**2 for i in range(len(BS_list)) ]
    # dYdz_1s_2 = [ BS_list[i].Gamma_break(z) * Y_1s[i] for i in range(len(BS_list)) ]
    # dYdz_1s_3 = [ BS_list[i].Gamma_ann_Hulthen_TA(z) * ( BS_list[i].Y_BS_eq(z) - Y_1s[i] ) for i in range(len(BS_list)) ]

    # dYdz_1s = prefactor_I * np.array( [ ( dYdz_1s_1[i] - dYdz_1s_2[i] + dYdz_1s_3[i] ) for i in range(len(BS_list)) ] )

    dYdz = np.concatenate( ([dYdz_DM], dYdz_1s) )

    return dYdz


def do_scan_Boltzmann_Hulthen_Network_1s(M_DM_scan, print_results):

    Omega_Boltzmann_Hulthen_Network_1s = []

    for M_DM in M_DM_scan:
        
        DM = Quintuplet_DM(M_DM)

        # # Variable order is: (Mass, gI, nE, l, I, nS, BS function, group factor constant)

        BS_1s1 = Quintuplet_BS(M_DM, gI_list[0], nE_list[0], l_list[0], Isp_list[0], nS_Quint, BS_func_list[0], gf_list[0])
        BS_1s3 = Quintuplet_BS(M_DM, gI_list[1], nE_list[1], l_list[1], Isp_list[1], nS_Quint, BS_func_list[1], gf_list[1])
        BS_1s5 = Quintuplet_BS(M_DM, gI_list[2], nE_list[2], l_list[2], Isp_list[2], nS_Quint, BS_func_list[2], gf_list[2])
        BS_list = [BS_1s1, BS_1s3, BS_1s5]

        bb = lambda z, Y: Boltzmann_Hulthen_Network_1s(DM, BS_list, z, Y)

        solution = NDE_solver(bb, len(BS_list) )
        
        omega_solution = Omega_DM(1.0, solution.y[0][-1], M_DM)

        Omega_Boltzmann_Hulthen_Network_1s.append([M_DM, omega_solution])

        if print_results == "y":
            print( [ M_DM, omega_solution ] )
        else: continue

    np.savetxt('./Results/Omega_Boltzmann_Hulthen_Network_1s.txt', Omega_Boltzmann_Hulthen_Network_1s)



############################################################################

# Define Boltzmann equation Network with all 1s, 2p, 2s BS

############################################################################ 


def Boltzmann_Hulthen_Network_1s2s2p(DM, BS_list, z, Y):


    Y_DM = Y[0]
    Y_1s1 = Y[1]
    Y_1s3 = Y[2]
    Y_1s5 = Y[3]
    Y_2s1 = Y[4]
    Y_2s3 = Y[5]
    Y_2s5 = Y[6]
    Y_2p1 = Y[7]
    Y_2p3 = Y[8]
    Y_2p5 = Y[9]
    Y_BS = [Y_1s1, Y_1s3, Y_1s5, Y_2s1, Y_2s3, Y_2s5, Y_2p1, Y_2p3, Y_2p5]

    ### Define DM mass ###

    M = DM.M

    ### Boltzmann equation for DM ###

    prefactor_DM = - DM.entropy(z)/( z * DM.Hubble(z) )

    dYdz_DM_1 = ( Y_DM**2 - DM.Yeq(z)**2 ) * DM.sv_production_Somm_Hulthen(z)
    dYdz_DM_2 = sigma0_prime(M) * np.sum( [ BS_list[i].bsf(z) * ( Y_DM**2 - Y_BS[i] * InverseYratio(BS_list[i], M, z) ) for i in range(len(BS_list)) ] )

    dYdz_DM = prefactor_DM * ( dYdz_DM_1 + dYdz_DM_2 )

    ### Boltzmann equation for BS ###

    prefactor_I = 1/( z * DM.Hubble(z) )

    dYdz_BS_1 = [ BS_list[i].Gamma_break_NoExp(z) * Yratio_NoExp(BS_list[i], M, z) * Y_DM**2 for i in range(len(BS_list)) ]
    dYdz_BS_2 = [ BS_list[i].Gamma_break(z) * Y_BS[i] for i in range(len(BS_list)) ]
    dYdz_BS_3 = [ BS_list[i].Gamma_ann_Hulthen_TA(z) * ( BS_list[i].Y_BS_eq(z) - Y_BS[i] ) for i in range(len(BS_list)) ]

    dYdz_BS = prefactor_I * np.array( [ ( dYdz_BS_1[i] - dYdz_BS_2[i] + dYdz_BS_3[i] ) for i in range(len(BS_list)) ] )

    dYdz = np.concatenate( ([dYdz_DM], dYdz_BS) )

    return dYdz



def do_scan_Boltzmann_Hulthen_Network_1s2s2p(M_DM_scan, print_results):

    Omega_Boltzmann_Hulthen_Network_1s2s2p = []

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

        bb = lambda z, Y: Boltzmann_Hulthen_Network_1s2s2p(DM, BS_list, z, Y)

        solution = NDE_solver(bb, len(BS_list) )
        
        try:
            omega_solution = Omega_DM(1.0, solution.y[0][-1], M_DM)
            Omega_Boltzmann_Hulthen_Network_1s2s2p.append([M_DM, omega_solution])

            if print_results == "y":
                print( [ M_DM, omega_solution ] )
            else: continue

        except IndexError:
            pass

    np.savetxt('./Results/Omega_Boltzmann_Hulthen_Network_1s2s2p.txt', Omega_Boltzmann_Hulthen_Network_1s2s2p)
