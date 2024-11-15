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


# ###############################
# # Effective couplings of BS
# ###############################

# def alpha_eff(Isp, n):
#     lambda_eff = (2*n*n - 1 - Isp*Isp)/8
#     return mpf(alpha*lambda_eff)


# ###############################

# # epsilon factor in Hulthen potential

# ############################### 


# def Hulthen_epsilon(M, z, Isp, nS):
#     return kappa * MWTSU2(M, z)/( alpha_eff(Isp, nS) * M )



# ###############################

# # Bindin energy of BS (for Coulomb potential)

# ###############################

# def binding_energy_BS(M, z, nE, l, Isp, nS):
#     prefactor = alpha_eff(Isp, nS)**2/(4 * nE**2)
#     H_epsilon = Hulthen_epsilon(M, z, Isp, nS)
#     corrections = 1 - nE**2 * H_epsilon - 0.53 * nE**2 * H_epsilon**2 * l * (l+1)
#     return float(prefactor * corrections**2)


# ###############################

# # BS number density

# ###############################


# def n_BS_eq(M, z, gI, nE, l, Isp, nS):
#     prefactor = gI * (2 * M**2/(2 * z * pi))**(3/2)
#     exponential = np.exp( -(2 - binding_energy_BS(M, z, nE, l, Isp, nS)) * z )
#     return mpf(prefactor * exponential)


# ###############################

# # BS abundance

# ###############################


# def Y_BS_eq(M, z, gI, nE, l, Isp, nS):
#     return n_BS_eq(M, z, gI, nE, l, Isp, nS)/entropy(M, z)


###############################

# Import thermal average of BSF cross sections !!! DIVIDED BY sigma0prime !!!

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


###############################

# Breaking rates of BS

###############################

def Gamma_break(M, z, gI, nE, l, Isp, nS, BSF_TA_xsec):
    pref_1 = gx**2/(2 * gI) * M**3 * sigma0_prime(M)
    pref_2 = (1/(z * 4 * pi))**(3/2)
    exponential = np.exp( - binding_energy_BS(M, z, nE, l, Isp, nS) * z )
    return pref_1 * pref_2 * exponential * BSF_TA_xsec(z)

###############################

# Annihilation rates of BS

###############################

def gamma_ann(group_factor, M):
    return group_factor * A25 * M

def Gamma_ann_Hulthen_TA(group_factor, M, z, Isp, nS):
    correction = 1 + Hulthen_epsilon(M, z, Isp, nS)**2
    return gamma_ann(group_factor, M) * correction * Kratio(z) 


######## Arrays of quanta = 9
# order: 1s1, 1s3, 1s5, 2s1, 2s3, 2s5, 2p1, 2p3, 2p5
# lambda_i, n ,l, g_I
lamba_eff = [6.0, 5.0, 3.0, 6.0, 5.0, 3.0, 6.0, 5.0, 3.0]
n_energy = [1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
l_angular = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
gI = [1.0, 9.0, 5.0, 1.0, 9.0, 5.0, 3.0, 3.0, 15.0]
## array of ann rate and transition rates

class Quintuplet_BS:
    """ Set of functions for the Bound States of Quintuplet DM """

    def __init__(self, DarkMatter_Mass, gI, nE, l, Isp, nS):
        self.M = DarkMatter_Mass
        self.gI = gI
        self.nE = nE
        self.l = l
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
    # Bindin energy of BS (for Coulomb potential)
    ###############################

    def binding_energy_BS(self, z):
        prefactor = self.alpha_eff()**2/(4 * self.nE**2)
        H_epsilon = self.Hulthen_epsilon(z)
        corrections = 1 - self.nE**2 * H_epsilon - 0.53 * self.nE**2 * H_epsilon**2 * self.l * (self.l+1)
        return float(prefactor * corrections**2)


    ###############################
    # BS number density
    ###############################


    def n_BS_eq(self, z):
        prefactor = self.gI * (2 * self.M**2/(2 * z * pi))**(3/2)
        exponential = np.exp( -(2 - self.binding_energy_BS(z)) * z )
        return mpf(prefactor * exponential)


    ###############################
    # BS abundance
    ###############################


    def Y_BS_eq(self, z):
        return self.n_BS_eq(z)/entropy(self.M, z)

    def func(self, b):
        return self.M * b


BS_1s1 = Quintuplet_BS(10 * TeV, gI[0], n_energy[0], l_angular[0], 1, 5)

print( BS_1s1.func(10) )
print( BS_1s1.alpha_eff() )
print( BS_1s1.binding_energy_BS(1e2) )
# Gann = [3240.0*A25*mdm, 15625.0/48.0*A25*mdm, 567.0/4.0*A25*mdm, 405.0*A25*mdm, 15625.0/384.0*A25*mdm, 567.0/32.0*A25*mdm, A25*a22*mdm, A25*a22*mdm, A25*a22*mdm ]
# Gdec = [0.0, 0.0, 0.0, s2tw*A25*s2tw*A2MZ*mdm, s2tw*A25*s2tw*A2MZ*mdm, s2tw*A25*s2tw*A2MZ*mdm, 2.0*s2tw*A25*mdm, 1.3*s2tw*A25*mdm, 0.2*s2tw*A25*mdm]
