import numpy as np
import sys
from mpmath import *
from scipy.integrate import solve_ivp
import time
from Basic_definitions import *
from Quintuplet_definitions import *
from BoundStates_definitions import *
from Free_Boltz_eqs import *
from Hulthen_eff_Boltz_eqs import *
from Hulthen_network_Boltz_eqs import *
mp.dps = 25; mp.pretty = True


############################################################################

# ASK FOR INPUTS

############################################################################

# print("Input number of masses to scan:")

n_scan = int(input("Input number of masses to scan: "))


M_DM_scan = np.linspace(3 * TeV, 15 * TeV, n_scan)


print("""

    Cases implemented:

    1 - Tree level (~0.6 s/mass)

    2 - Sommerfeld enhancement with Hulthen potential (~1.0 s/mass)

    3 - Effective Boltzmann equation with 1s bound states (~1.8 s/mass)

    4 - Effective Boltzmann equation with 1s, 2s, 2p bound states (~3.5 s/mass)

    5 - Network of Boltzmann equations with 1s1 bound states (~3.8 s/mass)

    6 - Network of Boltzmann equations with 1s bound states (~8.6 s/mass)

    7 - Network of Boltzmann equations with 1s, 2s, 2p bound states (~30 s/mass)

    0 - Solve all of them (could take some time)
    

    """)


which_case = int(input("Input case: "))

possible_cases = [0, 1, 2, 3, 4, 5, 6, 7]

if which_case not in possible_cases: 
    print("WRONG")
    quit()


print_results = input("Print results on terminal? y, n (default): ")
if len(print_results) == 0: print_results = "n"



# ############################################################################

# # Do selected scan

# ############################################################################ 

############################################################################

start = time.time()

############################################################################


if which_case == 0:
    # Scan with free Boltzmann equation
    do_scan_Boltzmann_Tree(M_DM_scan, print_results)
    # Scan with free Boltzmann equation, using Hulthen potential
    do_scan_Boltzmann_Hulthen_free(M_DM_scan, print_results)
    # Scan with effective Boltzmann equation, using all 1s BS
    do_scan_Boltzmann_Hulthen_1s_effective(M_DM_scan, print_results)
    # Scan with effective Boltzmann equation, using all 1s, 2s, 2p BS
    do_scan_Boltzmann_Hulthen_1s2s2p_effective(M_DM_scan, print_results)
    # Scan with network of Boltzmann equations, using 1s1 BS
    do_scan_Boltzmann_Hulthen_Network_1s1(M_DM_scan, print_results)
    # Scan with network of Boltzmann equations, using 1s BS
    do_scan_Boltzmann_Hulthen_Network_1s(M_DM_scan, print_results)
    # Scan with network of Boltzmann equations, using all 1s, 2s, 2p BS
    do_scan_Boltzmann_Hulthen_Network_1s2s2p(M_DM_scan, print_results)
elif which_case == 1:
    do_scan_Boltzmann_Tree(M_DM_scan, print_results)
elif which_case == 2:
    do_scan_Boltzmann_Hulthen_free(M_DM_scan, print_results)
elif which_case == 3:
    do_scan_Boltzmann_Hulthen_1s_effective(M_DM_scan, print_results)
elif which_case == 4:
    do_scan_Boltzmann_Hulthen_1s2s2p_effective(M_DM_scan, print_results)
elif which_case == 5:
    do_scan_Boltzmann_Hulthen_Network_1s1(M_DM_scan, print_results)
elif which_case == 6:
    do_scan_Boltzmann_Hulthen_Network_1s(M_DM_scan, print_results)
elif which_case == 7:
    do_scan_Boltzmann_Hulthen_Network_1s2s2p(M_DM_scan, print_results)


############################################################################

end = time.time()
print("Total computing time: ", end - start)

############################################################################