import numpy as np
from mpmath import *
from scipy import integrate, interpolate
import time
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

# Import various quantities computed numerically

############################################################################ 

###############################
# Bessel functions K1 and K2
###############################

K1_hiprec = np.loadtxt('Num_lists/K1_hiprec.txt')
K2_hiprec = np.loadtxt('Num_lists/K2_hiprec.txt')

k1_interp = interpolate.interp1d( K1_hiprec[:, 0], K1_hiprec[:, 1] )
k2_interp = interpolate.interp1d( K2_hiprec[:, 0], K2_hiprec[:, 1] )

def K1(z):
    return 10**mpf( float( k1_interp(z) ) )

def K2(z):
    return 10**mpf( float( k2_interp(z) ) )

###############################
# K2(z)/z (for number density)
###############################

nDM_factor = np.loadtxt('Num_lists/nDMeq_z.txt')

nDM_factor_interp = interpolate.interp1d( nDM_factor[:, 0], nDM_factor[:, 1] )

def K2_over_z(z):
    return 10**mpf( float( nDM_factor_interp(z) ) )


###############################
# K2(z)*z^2 (for equilibrium abundance)
###############################

Y_factor = np.loadtxt('Num_lists/YDMeq_z.txt')

Y_factor_interp = interpolate.interp1d( Y_factor[:, 0], Y_factor[:, 1] )

def K2_times_z2(z):
    return 10**mpf( float( Y_factor_interp(z) ) )


###############################
# Thermal average of BSF cross sections !!! DIVIDED BY sigma0prime !!!
###############################

BSF_1s1_TA_num = np.loadtxt('Num_lists/BSF_1s1_TA.txt')
BSF_1s3_TA_num = np.loadtxt('Num_lists/BSF_1s3_TA.txt')
BSF_1s5_TA_num = np.loadtxt('Num_lists/BSF_1s5_TA.txt')
BSF_2p1_TA_num = np.loadtxt('Num_lists/BSF_2p1_TA.txt')
BSF_2p3_TA_num = np.loadtxt('Num_lists/BSF_2p3_TA.txt')
BSF_2p5_TA_num = np.loadtxt('Num_lists/BSF_2p5_TA.txt')
BSF_2s1_TA_num = np.loadtxt('Num_lists/BSF_2s1_TA.txt')
BSF_2s3_TA_num = np.loadtxt('Num_lists/BSF_2s3_TA.txt')
BSF_2s5_TA_num = np.loadtxt('Num_lists/BSF_2s5_TA.txt')

BSF_1s1_TA_num_interp = interpolate.interp1d( BSF_1s1_TA_num[:, 0], BSF_1s1_TA_num[:, 1] )
BSF_1s3_TA_num_interp = interpolate.interp1d( BSF_1s3_TA_num[:, 0], BSF_1s3_TA_num[:, 1] )
BSF_1s5_TA_num_interp = interpolate.interp1d( BSF_1s5_TA_num[:, 0], BSF_1s5_TA_num[:, 1] )
BSF_2p1_TA_num_interp = interpolate.interp1d( BSF_2p1_TA_num[:, 0], BSF_2p1_TA_num[:, 1] )
BSF_2p3_TA_num_interp = interpolate.interp1d( BSF_2p3_TA_num[:, 0], BSF_2p3_TA_num[:, 1] )
BSF_2p5_TA_num_interp = interpolate.interp1d( BSF_2p5_TA_num[:, 0], BSF_2p5_TA_num[:, 1] )
BSF_2s1_TA_num_interp = interpolate.interp1d( BSF_2s1_TA_num[:, 0], BSF_2s1_TA_num[:, 1] )
BSF_2s3_TA_num_interp = interpolate.interp1d( BSF_2s3_TA_num[:, 0], BSF_2s3_TA_num[:, 1] )
BSF_2s5_TA_num_interp = interpolate.interp1d( BSF_2s5_TA_num[:, 0], BSF_2s5_TA_num[:, 1] )

def BSF_1s1_TA(z):
    return 10**mpf( float( BSF_1s1_TA_num_interp(z) ) )

def BSF_1s3_TA(z):
    return 10**mpf( float( BSF_1s3_TA_num_interp(z) ) )

def BSF_1s5_TA(z):
    return 10**mpf( float( BSF_1s5_TA_num_interp(z) ) )

def BSF_2p1_TA(z):
    return 10**mpf( float( BSF_2p1_TA_num_interp(z) ) )

def BSF_2p3_TA(z):
    return 10**mpf( float( BSF_2p3_TA_num_interp(z) ) )

def BSF_2p5_TA(z):
    return 10**mpf( float( BSF_2p5_TA_num_interp(z) ) )

def BSF_2s1_TA(z):
    return 10**mpf( float( BSF_2s1_TA_num_interp(z) ) )

def BSF_2s3_TA(z):
    return 10**mpf( float( BSF_2s3_TA_num_interp(z) ) )

def BSF_2s5_TA(z):
    return 10**mpf( float( BSF_2s5_TA_num_interp(z) ) )





