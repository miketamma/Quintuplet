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

# ASK FOR INPUTS

############################################################################

# print("Input number of masses to scan:")

n_scan = int(input("Input number of masses to scan: "))


M_DM_scan = np.linspace(3 * TeV, 14 * TeV, n_scan)


print("""

    Cases implemented:

    1 - Tree level (~0.6 s/mass)

    2 - Sommerfeld enhancement with Hulthen potential (~1.0 s/mass)

    3 - Effective Boltzmann equation with 1s bound states (~1.8 s/mass)

    4 - Effective Boltzmann equation with 1s, 2s, 2p bound states (~3.5 s/mass)

    5 - Network of Boltzmann equations with 1s1 bound states (~3.8 s/mass)

    6 - Network of Boltzmann equations with 1s bound states (~8.6 s/mass)

    0 - Solve all of them (could take some time)
    

    """)


which_case = int(input("Input case: "))

possible_cases = [0, 1, 2, 3, 4, 5, 6]

if which_case not in possible_cases: 
    print("WRONG")
    quit()


print_results = input("Print results on terminal? y, n (default): ")
if len(print_results) == 0: print_results = "n"




############################################################################

start = time.time()

############################################################################



############################################################################

# Numerical ODE solver, n is the number of variables in the ODE

############################################################################

def NDE_solver(Boltz, n):
    y0 = np.full( n, [1.e-5,])
    teval = np.logspace(np.log10(3.5), 5.9, 100)
    tspan=[3.0, 1.e6]
    
    return solve_ivp(Boltz, tspan, y0 = y0, method = 'Radau', t_eval = teval, atol=1.e-20)


############################################################################

# Define Boltzmann equation without BS (Tree level)

############################################################################ 

def Boltzmann_Tree(DM, z, Y):

    prefactor = - DM.entropy(z)/( z * DM.Hubble(z) )

    abundances = Y**2 - DM.Yeq(z)**2

    return prefactor * DM.sv_production_NoSomm(z) * abundances


def do_scan_Boltzmann_Tree():

    Omega_Boltzmann_Tree = []

    for M_DM in M_DM_scan:
        
        DM = Quintuplet_DM(M_DM)

        bb = lambda z, Y: Boltzmann_Tree(DM, z, Y)

        solution = NDE_solver(bb, 1)

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


def do_scan_Boltzmann_Hulthen_free():

    Omega_Boltzmann_Hulthen_free = []

    for M_DM in M_DM_scan:
        
        DM = Quintuplet_DM(M_DM)

        bb = lambda z, Y: Boltzmann_Hulthen_free(DM, z, Y)

        solution = NDE_solver(bb, 1)

        omega_solution = Omega_DM(1.0, solution.y[0][-1], M_DM)

        Omega_Boltzmann_Hulthen_free.append([M_DM, omega_solution])

        if print_results == "y":
            print( [ M_DM, omega_solution ] )
        else: continue

    np.savetxt('./Results/Omega_Boltzmann_Hulthen_free.txt', Omega_Boltzmann_Hulthen_free)



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

def do_scan_Boltzmann_Hulthen_1s_effective():

    Omega_Boltzmann_Hulthen_1s_effective = []

    for M_DM in M_DM_scan:
        
        DM = Quintuplet_DM(M_DM)

        # # Variable order is: (Mass, gI, nE, l, I, nS, BS function, group factor constant)

        BS_1s1 = Quintuplet_BS(M_DM, gI_list[0], nE_list[0], l_list[0], Isp_list[0], nS_Quint, BS_func_list[0], gf_list[0])
        BS_1s3 = Quintuplet_BS(M_DM, gI_list[1], nE_list[1], l_list[1], Isp_list[1], nS_Quint, BS_func_list[1], gf_list[1])
        BS_1s5 = Quintuplet_BS(M_DM, gI_list[2], nE_list[2], l_list[2], Isp_list[2], nS_Quint, BS_func_list[2], gf_list[2])

        BS_1s_list = [BS_1s1, BS_1s3, BS_1s5]

        bb = lambda z, Y: Boltzmann_Hulthen_effective(DM, BS_1s_list, z, Y)

        solution = NDE_solver(bb, 1)
        
        omega_solution = Omega_DM(1.0, solution.y[0][-1], M_DM)

        Omega_Boltzmann_Hulthen_1s_effective.append([M_DM, omega_solution])

        if print_results == "y":
            print( [ M_DM, omega_solution ] )
        else: continue

    np.savetxt('./Results/Omega_Boltzmann_Hulthen_1s_effective.txt', Omega_Boltzmann_Hulthen_1s_effective)



############################################################################

# Define Boltzmann equation with 1s, 2s, 2p BS (effective equation)

############################################################################ 


# BS_1s1 = Quintuplet_BS(10 * TeV, gI_list[0], nE_list[0], l_list[0], Isp_list[0], nS_Quint, BS_func_list[0], gf_list[0])
# BS_1s3 = Quintuplet_BS(10 * TeV, gI_list[0], nE_list[0], l_list[0], Isp_list[0], nS_Quint, BS_func_list[0], gf_list[0])
# BS_1s5 = Quintuplet_BS(10 * TeV, gI_list[0], nE_list[0], l_list[0], Isp_list[0], nS_Quint, BS_func_list[0], gf_list[0])

# BS_1s_list = [BS_1s1,BS_1s3,BS_1s5] 

# print( BS_eff(BS_1s1, 10 * TeV, 1e3) )
# print( np.sum([BS_eff(BS, 10 * TeV, 1e3) for BS in BS_1s_list] ) )



def do_scan_Boltzmann_Hulthen_1s2s2p_effective():

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

        solution = NDE_solver(bb, 1)
        
        omega_solution = Omega_DM(1.0, solution.y[0][-1], M_DM)

        Omega_Boltzmann_Hulthen_1s2s2p_effective.append([M_DM, omega_solution])

        if print_results == "y":
            print( [ M_DM, omega_solution ] )
        else: continue

    np.savetxt('./Results/Omega_Boltzmann_Hulthen_1s2s2p_effective.txt', Omega_Boltzmann_Hulthen_1s2s2p_effective)



############################################################################

# Define ratios for the network cases

############################################################################

def Yratio(BS, M, z):

    prefactor = ( pi**(5/2) ) * 32/45

    prefactor_g = BS.gI * gs_star(M/z)/( gx**2 )

    factor_z = ( np.exp( z * BS.binding_energy_BS(z) ) )/( z**(3/2) )

    return prefactor * prefactor_g * factor_z

def Yratio_NoExp(BS, M, z):

    prefactor = ( pi**(5/2) ) * 32/45

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
    dYdz_DM_2 = sigma0_prime(M) * BS.bsf(z) * ( Y_DM**2 - Y_1s1/Yratio(BS, M, z) )

    dYdz_DM = prefactor_DM * ( dYdz_DM_1 + dYdz_DM_2 )

    ### Boltzmann equation for BS ###

    prefactor_I = 1/( z * DM.Hubble(z) )

    dYdz_1s1_1 = BS.Gamma_break_NoExp(z) * Yratio_NoExp(BS, M, z) * Y_DM**2
    dYdz_1s1_2 = - BS.Gamma_break(z) * Y_1s1
    dYdz_1s1_3 = BS.Gamma_ann_Hulthen_TA(z) * ( BS.Y_BS_eq(z) - Y_1s1 )

    dYdz_1s1 = prefactor_I * ( dYdz_1s1_1 + dYdz_1s1_2 + dYdz_1s1_3)

    return [dYdz_DM, dYdz_1s1]


def do_scan_Boltzmann_Hulthen_Network_1s1():

    Omega_Boltzmann_Hulthen_Network_1s1 = []

    for M_DM in M_DM_scan:
        
        DM = Quintuplet_DM(M_DM)

        # # Variable order is: (Mass, gI, nE, l, I, nS, BS function, group factor constant)

        BS_1s1 = Quintuplet_BS(M_DM, gI_list[0], nE_list[0], l_list[0], Isp_list[0], nS_Quint, BS_func_list[0], gf_list[0])

        bb = lambda z, Y: Boltzmann_Hulthen_Network_1s1(DM, BS_1s1, z, Y)

        solution = NDE_solver(bb, 2)
        
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
    prefactor_I = 1/( z * DM.Hubble(z) )

    dYdz_DM_1 = ( Y_DM**2 - DM.Yeq(z)**2 ) * DM.sv_production_Somm_Hulthen(z)
    dYdz_DM_2 = sigma0_prime(M) * np.sum( [ BS_list[i].bsf(z) * ( Y_DM**2 - Y_1s[i]/Yratio(BS_list[i], M, z) ) for i in range(len(BS_list)) ] )

    dYdz_DM = prefactor_DM * ( dYdz_DM_1 + dYdz_DM_2 )

    ### Boltzmann equation for BS ###

    prefactor_I = 1/( z * DM.Hubble(z) )

    dYdz_1s_1 = [ BS.Gamma_break_NoExp(z) * Yratio_NoExp(BS, M, z) * Y_DM**2 for BS in BS_list ]
    dYdz_1s_2 = [- BS_list[i].Gamma_break(z) * Y_1s[i] for i in range(len(BS_list)) ]
    dYdz_1s_3 = [ BS_list[i].Gamma_ann_Hulthen_TA(z) * ( BS_list[i].Y_BS_eq(z) - Y_1s[i] ) for i in range(len(BS_list)) ]

    dYdz_1s = prefactor_I * np.array( [ ( dYdz_1s_1[i] + dYdz_1s_2[i] + dYdz_1s_3[i] ) for i in range(len(BS_list)) ] )

    dYdz = np.concatenate( ([dYdz_DM], dYdz_1s) )

    return dYdz

# DM = Quintuplet_DM(10*TeV)
# BS_1s1 = Quintuplet_BS(10*TeV, gI_list[0], nE_list[0], l_list[0], Isp_list[0], nS_Quint, BS_func_list[0], gf_list[0])
# BS_1s3 = Quintuplet_BS(10*TeV, gI_list[1], nE_list[1], l_list[1], Isp_list[1], nS_Quint, BS_func_list[1], gf_list[1])
# BS_1s5 = Quintuplet_BS(10*TeV, gI_list[2], nE_list[2], l_list[2], Isp_list[2], nS_Quint, BS_func_list[2], gf_list[2])
# BS_list = [BS_1s1, BS_1s3, BS_1s5]

# print( Boltzmann_Hulthen_Network_1s(DM, BS_list, 100, [0,0,0,0]) )

def do_scan_Boltzmann_Hulthen_Network_1s():

    Omega_Boltzmann_Hulthen_Network_1s = []

    for M_DM in M_DM_scan:
        
        DM = Quintuplet_DM(M_DM)

        # # Variable order is: (Mass, gI, nE, l, I, nS, BS function, group factor constant)

        BS_1s1 = Quintuplet_BS(10*TeV, gI_list[0], nE_list[0], l_list[0], Isp_list[0], nS_Quint, BS_func_list[0], gf_list[0])
        BS_1s3 = Quintuplet_BS(10*TeV, gI_list[1], nE_list[1], l_list[1], Isp_list[1], nS_Quint, BS_func_list[1], gf_list[1])
        BS_1s5 = Quintuplet_BS(10*TeV, gI_list[2], nE_list[2], l_list[2], Isp_list[2], nS_Quint, BS_func_list[2], gf_list[2])
        BS_list = [BS_1s1, BS_1s3, BS_1s5]

        bb = lambda z, Y: Boltzmann_Hulthen_Network_1s(DM, BS_list, z, Y)

        solution = NDE_solver(bb, 4)
        
        omega_solution = Omega_DM(1.0, solution.y[0][-1], M_DM)

        Omega_Boltzmann_Hulthen_Network_1s.append([M_DM, omega_solution])

        if print_results == "y":
            print( [ M_DM, omega_solution ] )
        else: continue

    np.savetxt('./Results/Omega_Boltzmann_Hulthen_Network_1s.txt', Omega_Boltzmann_Hulthen_Network_1s)


############################################################################

# Do selected scan

############################################################################ 

if which_case == 0:
    # Scan with free Boltzmann equation
    do_scan_Boltzmann_Tree()
    # Scan with free Boltzmann equation, using Hulthen potential
    do_scan_Boltzmann_Hulthen_free()
    # Scan with effective Boltzmann equation, using all 1s BS
    do_scan_Boltzmann_Hulthen_1s_effective()
    # Scan with effective Boltzmann equation, using all 1s, 2s, 2p BS
    do_scan_Boltzmann_Hulthen_1s2s2p_effective()
    # Scan with network of Boltzmann equations, using 1s1 BS
    do_scan_Boltzmann_Hulthen_Network_1s1()
    # Scan with network of Boltzmann equations, using 1s BS
    do_scan_Boltzmann_Hulthen_Network_1s()
elif which_case == 1:
    do_scan_Boltzmann_Tree()
elif which_case == 2:
    do_scan_Boltzmann_Hulthen_free()
elif which_case == 3:
    do_scan_Boltzmann_Hulthen_1s_effective()
elif which_case == 4:
    do_scan_Boltzmann_Hulthen_1s2s2p_effective()
elif which_case == 5:
    do_scan_Boltzmann_Hulthen_Network_1s1()
elif which_case == 6:
    do_scan_Boltzmann_Hulthen_Network_1s()



end = time.time()
print("Total computing time: ", end - start)