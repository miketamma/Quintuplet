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
  return re( sqrt(1 - ((M/z)**2)/(Tc**2) ) )

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

# print( SommHulthen_S1(kappa * MWTSU2(10 * TeV, 1e6)/( alpha * lam[0] * 10 * TeV ), 1e6) )

class Quintuplet_DM:
    """ Set of functions for the Quintuplet DM """

    def __init__(self, DarkMatter_Mass):
        self.M = DarkMatter_Mass


    ###############################
    # Effective couplings for Quintuplet
    ###############################

    def lam(self):
        return [6,5,3]

    ###############################
    # epsilon factor in Hulther potential
    ###############################

    def h_eps(self, z, i):
        return kappa * MWTSU2(self.M, z)/( alpha * self.lam()[i] * self.M )

    ###############################
    # entropy and Hubble for DM
    ###############################

    def entropy(self, z):
        return 2*pi**2/45 * ( sqrtg_interp(self.M/z)**2 ) * (self.M/z)**3

    def Hubble(self, z):
        return 2 * re( sqrt( pi**3/45 * ( sqrtg_interp(self.M/z)**2 ) ) ) * self.M**2/(Mpl * z**2)

    ###############################
    # Equilibrium abundance for DM
    ###############################

    def Yeq(self, z):
        pref = gx * 45 / ( 2 * pi**2 * (2 * pi)**(3/2) * sqrtg_interp(self.M/z)**2 )
        return pref * K2_times_z2(z)

    ###############################
    # Production cross section with Sommerfeld enhancement (Hulthen potential)
    ###############################

    def sv_production(self, z):
        S1 = 16/69 * SommHulthen_S1(self.h_eps(z, 0), z)
        S3 = 25/69 * SommHulthen_S3(self.h_eps(z, 1), z)
        S5 = 28/69 * SommHulthen_S5(self.h_eps(z, 2), z)
        return sigma0(self.M) * (S1 + S3 + S5)

