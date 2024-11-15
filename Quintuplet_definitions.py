import numpy as np
from mpmath import *
from scipy import integrate, interpolate
import time
from Basic_definitions import *
mp.dps = 25; mp.pretty = True
######################################


############################################################################

# W mass as function of temperature

############################################################################ 

def Rvev(M, z):
  return np.real( np.sqrt(1 - ((M/z)**2)/(Tc**2) ) )

def MWSU2(M, z):
  return np.sqrt( 11/6 * ((M/z)**2) * alpha * 4 * pi )

def MWTSU2(M, z):
  return np.sqrt( MW**2 * Rvev(M,z)**2 + MWSU2(M,z)**2 )


############################################################################

# s-wave annihilation cross section for Quintuplet

############################################################################ 

def sigma0(M):
  return 207.0/20.0*pi*A2MZ**2.0/M**2.0 

def sigma0_prime(M):
  return pi*A2MZ**2.0/M**2.0 


############################################################################

# Definitions for BS quantities

############################################################################


###############################

# Import thermal average of Sommerfeld enhancement !!! DIVIDED BY sigma0 !!!

###############################


epsilon, zeta, SommHulthen_S1_num, SommHulthen_S3_num, SommHulthen_S5_num = np.loadtxt('Num_lists/SommHueltenTA.txt', unpack = True)


SommHulthen_S1_interp = interpolate.NearestNDInterpolator(list(zip(epsilon, zeta)), SommHulthen_S1_num)
SommHulthen_S3_interp = interpolate.NearestNDInterpolator(list(zip(epsilon, zeta)), SommHulthen_S3_num)
SommHulthen_S5_interp = interpolate.NearestNDInterpolator(list(zip(epsilon, zeta)), SommHulthen_S5_num)

def SommHulthen_S1(eps, z):
    return 10**mpf( float( SommHulthen_S1_interp(eps, z) ) )

def SommHulthen_S3(eps, z):
    return 10**mpf( float( SommHulthen_S3_interp(eps, z) ) )

def SommHulthen_S5(eps, z):
    return 10**mpf( float( SommHulthen_S5_interp(eps, z) ) )



class Quintuplet_DM:
    """ Set of functions for the Quintuplet DM """

    def __init__(self, DarkMatter_Mass, Isp, nS):
        self.M = DarkMatter_Mass
        self.Isp = Isp
        self.nS = nS

    ###############################
    # Effective couplings of BS
    ###############################

    def alpha_eff(self):
        lambda_eff = (2 * self.nS**2 - 1 - self.Isp**2)/8
        return mpf(alpha*lambda_eff)


    ###############################
    # epsilon factor in Hulthen potential
    ############################### 

    def Hulthen_epsilon(self, z):
        return kappa * MWTSU2(self.M, z)/( self.alpha_eff() * self.M )

    ###############################
    # Production cross section with Sommerfeld enhancement (Hulthen potential)
    ###############################

    def sv_production(self, z):
        S1 = 16/69 * SommHulthen_S1(self.Hulthen_epsilon(z), z)
        S3 = 25/69 * SommHulthen_S3(self.Hulthen_epsilon(z), z)
        S5 = 28/69 * SommHulthen_S5(self.Hulthen_epsilon(z), z)
        return sigma0(self.M) * (S1 + S3 + S5)



# BS_1s1 = Quintuplet_BS(10 * TeV, gI[0], n_energy[0], l_angular[0], 1, 5, BSF_1s1_TA, 3240.0)

# print( BS_1s1.alpha_eff() )
# print( BS_1s1.n_BS_eq(1e6) )
# print( BS_1s1.gamma_ann() )