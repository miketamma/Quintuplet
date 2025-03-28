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
import matplotlib.pyplot as plt
mp.dps = 25; mp.pretty = True


M_DM = 10 * TeV


DM = Quintuplet_DM(M_DM)

# # Variable order is: (Mass, gI, nE, l, I, nS, BS function, group factor constant)

BS_1s1 = Quintuplet_BS(M_DM, gI_list[0], nE_list[0], l_list[0], Isp_list[0], nS_Quint, BS_func_list[0], gf_list[0])
BS_1s3 = Quintuplet_BS(M_DM, gI_list[1], nE_list[1], l_list[1], Isp_list[1], nS_Quint, BS_func_list[1], gf_list[1])
BS_1s5 = Quintuplet_BS(M_DM, gI_list[2], nE_list[2], l_list[2], Isp_list[2], nS_Quint, BS_func_list[2], gf_list[2])

BS_1s_list = [BS_1s1, BS_1s3, BS_1s5]

bb = lambda z, Y: Boltzmann_Hulthen_effective(DM, BS_1s_list, z, Y)
solution = NDE_solver(bb, 0)

bb1 = lambda z, Y: Boltzmann_Hulthen_Network_1s(DM, BS_1s_list, z, Y)
solution1 = NDE_solver(bb1, len(BS_1s_list) )

eff = [[solution.t[i], solution.y[0][i]] for i in range(len(solution.y[0]))]

eff1 = [[solution1.t[i], solution1.y[0][i]] for i in range(len(solution1.y[0]))]
eff1_1s1 = [[solution1.t[i], solution1.y[1][i]] for i in range(len(solution1.y[1]))]
eff1_1s3 = [[solution1.t[i], solution1.y[2][i]] for i in range(len(solution1.y[2]))]
eff1_1s5 = [[solution1.t[i], solution1.y[3][i]] for i in range(len(solution1.y[3]))]

# np.savetxt('./Results/db_eff.txt', eff)
# np.savetxt('./Results/db_eff1.txt', eff1)
# np.savetxt('./Results/db_eff1_1s1.txt', eff1_1s1)
# np.savetxt('./Results/db_eff1_1s3.txt', eff1_1s3)
# np.savetxt('./Results/db_eff1_1s5.txt', eff1_1s5)


# plt.subplot(2,1,1)
plt.plot([ eff[i][0] for i in range(len(eff)) ], [ eff[i][1] for i in range(len(eff)) ], color='blue', lw=2)
plt.plot([ eff1[i][0] for i in range(len(eff1)) ], [ eff1[i][1] for i in range(len(eff1)) ], color='red', lw=2)
# plt.plot([ eff1_1s1[i][0] for i in range(len(eff1_1s1)) ], [ eff1_1s1[i][1] for i in range(len(eff1_1s1)) ], '-.', color='orange', lw=2)
# plt.plot([ eff1_1s3[i][0] for i in range(len(eff1_1s3)) ], [ eff1_1s3[i][1] for i in range(len(eff1_1s3)) ], '.', color='purple', lw=2)
# plt.plot([ eff1_1s5[i][0] for i in range(len(eff1_1s5)) ], [ eff1_1s5[i][1] for i in range(len(eff1_1s5)) ], '-', color='green', lw=2)
plt.yscale('log')
plt.xscale('log')
plt.show()

