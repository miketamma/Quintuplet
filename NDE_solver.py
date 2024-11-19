import numpy as np
import sys
from mpmath import *
from scipy.integrate import solve_ivp
import time
from Basic_definitions import *
from Quintuplet_definitions import *
from BoundStates_definitions import *
mp.dps = 25; mp.pretty = True
######################################

start = time.time()

############################################################################

# DM mass range in TeV

############################################################################

M_DM_scan = np.linspace(3 * TeV, 14 * TeV, 100)


# DM = Quintuplet_DM(M_DM)

# # Variable order is: (Mass, gI, nE, l, I, nS, BS function, group factor constant)

# BS_1s1 = Quintuplet_BS(M_DM, gI_list[0], nE_list[0], l_list[0], Isp_list[0], nS_Quint, BS_func_list[0], gf_list[0])
# BS_1s3 = Quintuplet_BS(M_DM, gI_list[1], nE_list[1], l_list[1], Isp_list[1], nS_Quint, BS_func_list[1], gf_list[1])
# BS_1s5 = Quintuplet_BS(M_DM, gI_list[2], nE_list[2], l_list[2], Isp_list[2], nS_Quint, BS_func_list[2], gf_list[2])

# BS_2s1 = Quintuplet_BS(M_DM, gI_list[3], nE_list[3], l_list[3], Isp_list[3], nS_Quint, BS_func_list[3], gf_list[3])
# BS_2s3 = Quintuplet_BS(M_DM, gI_list[4], nE_list[4], l_list[4], Isp_list[4], nS_Quint, BS_func_list[4], gf_list[4])
# BS_2s5 = Quintuplet_BS(M_DM, gI_list[5], nE_list[5], l_list[5], Isp_list[5], nS_Quint, BS_func_list[5], gf_list[5])

# BS_2p1 = Quintuplet_BS(M_DM, gI_list[6], nE_list[6], l_list[6], Isp_list[6], nS_Quint, BS_func_list[6], gf_list[6])
# BS_2p3 = Quintuplet_BS(M_DM, gI_list[7], nE_list[7], l_list[7], Isp_list[7], nS_Quint, BS_func_list[7], gf_list[7])
# BS_2p5 = Quintuplet_BS(M_DM, gI_list[8], nE_list[8], l_list[8], Isp_list[8], nS_Quint, BS_func_list[8], gf_list[8])


# BS_1s_list = [BS_1s1, BS_1s3, BS_1s5]
# BS_2s_list = [BS_2s1, BS_2s3, BS_2s5]
# BS_2p_list = [BS_2p1, BS_2p3, BS_2p5]

############################################################################

# Numerical ODE solver, n is the number of variables in the ODE

############################################################################

def NDE_solver(Boltz, n):
    y0 = np.full( (n), [1.e-30,])
    teval = np.logspace(np.log10(3.5), 5.9, 100)
    tspan=[3.0, 1.e6]
    
    return solve_ivp(Boltz, tspan, y0 = y0, method = 'Radau', t_eval = teval, atol=1.e-20)


############################################################################

# Define Boltzmann equation without BS

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

        print( Omega_DM(1.0, solution.y[0][-1], M_DM) )
        Omega_Boltzmann_Tree.append([M_DM, omega_solution])

    np.savetxt('./Results/Omega_Boltzmann_Tree.txt', Omega_Boltzmann_Tree)



############################################################################

# Define Boltzmann equation without BS

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

        solution_free = NDE_solver(bb, 1)

        omega_solution = Omega_DM(1.0, solution_free.y[0][-1], M_DM)

        print( Omega_DM(1.0, solution_free.y[0][-1], M_DM) )
        Omega_Boltzmann_Hulthen_free.append([M_DM, omega_solution])

    np.savetxt('./Results/Omega_Boltzmann_Hulthen_free.txt', Omega_Boltzmann_Hulthen_free)



############################################################################

# Define Boltzmann equation with 1s BS (effective equation)

############################################################################ 


def BS_eff(BS, M, z):

    first_term = 1/BS.bsf(z)
    second_term_1 = gx**2 * sigma0_prime(M) * M**2/( 2 * BS.gI * BS.gamma_ann() )
    second_term_2 = ( 1 / (4*pi*z) )**(3/2)
    second_term_3 = np.exp( - z * BS.binding_energy_BS(z) )

    return sigma0_prime(M) * 1/(first_term + second_term_1 * second_term_2 * second_term_3)

    # return [sigma0_prime(M)/( * np.exp( - z * BS_list[i].binding_energy_BS(z)) ) for i in range(len(BS_list))]

# BS_1s1 = Quintuplet_BS(10 * TeV, gI_list[0], nE_list[0], l_list[0], Isp_list[0], nS_Quint, BS_func_list[0], gf_list[0])
# BS_1s3 = Quintuplet_BS(10 * TeV, gI_list[0], nE_list[0], l_list[0], Isp_list[0], nS_Quint, BS_func_list[0], gf_list[0])
# BS_1s5 = Quintuplet_BS(10 * TeV, gI_list[0], nE_list[0], l_list[0], Isp_list[0], nS_Quint, BS_func_list[0], gf_list[0])

# BS_1s_list = [BS_1s1,BS_1s3,BS_1s5] 

# print( BS_eff(BS_1s1, 10 * TeV, 1e3) )
# print( np.sum([BS_eff(BS, 10 * TeV, 1e3) for BS in BS_1s_list] ) )

def Boltzmann_Hulthen_1s_effective(DM, BS_list, z, Y):

    M = DM.M

    prefactor = - DM.entropy(z)/( z * DM.Hubble(z) )

    abundances = Y**2 - DM.Yeq(z)**2

    sv_eff = DM.sv_production_Somm_Hulthen(z) + np.sum([BS_eff(BS, M, z) for BS in BS_list])

    return prefactor * sv_eff * abundances



def do_scan_Boltzmann_Hulthen_1s_effective():

    Omega_Boltzmann_Hulthen_1s_effective = []

    for M_DM in M_DM_scan:
        
        DM = Quintuplet_DM(M_DM)

        # # Variable order is: (Mass, gI, nE, l, I, nS, BS function, group factor constant)

        BS_1s1 = Quintuplet_BS(M_DM, gI_list[0], nE_list[0], l_list[0], Isp_list[0], nS_Quint, BS_func_list[0], gf_list[0])
        BS_1s3 = Quintuplet_BS(M_DM, gI_list[1], nE_list[1], l_list[1], Isp_list[1], nS_Quint, BS_func_list[1], gf_list[1])
        BS_1s5 = Quintuplet_BS(M_DM, gI_list[2], nE_list[2], l_list[2], Isp_list[2], nS_Quint, BS_func_list[2], gf_list[2])

        BS_1s_list = [BS_1s1, BS_1s3, BS_1s5]

        bb = lambda z, Y: Boltzmann_Hulthen_1s_effective(DM, BS_1s_list, z, Y)

        solution = NDE_solver(bb, 1)
        
        omega_solution = Omega_DM(1.0, solution.y[0][-1], M_DM)

        print( Omega_DM(1.0, solution.y[0][-1], M_DM) )
        Omega_Boltzmann_Hulthen_1s_effective.append([M_DM, omega_solution])

    np.savetxt('./Results/Omega_Boltzmann_Hulthen_1s_effective.txt', Omega_Boltzmann_Hulthen_1s_effective)


############################################################################

# Do selected scan

############################################################################ 

# Scan with free Boltzmann equation
do_scan_Boltzmann_Tree()

# Scan with free Boltzmann equation, using Hulthen potential
# do_scan_Boltzmann_Hulthen_free()


# Scan with effective Boltzmann equation, using all 1s BS
# do_scan_Boltzmann_Hulthen_1s_effective()

end = time.time()
print(end - start)