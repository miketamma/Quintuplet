import numpy as np
from mpmath import *
from scipy import integrate, interpolate
import time
from Basic_definitions import *
from Quintuplet_definitions import *
mp.dps = 25; mp.pretty = True
######################################


############################################################################

# Definitions for BS quantities

############################################################################


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
# Effective couplings of BS
###############################

def alpha_eff(Isp, nS):
    lambda_eff = (2 * nS**2 - 1 - Isp**2)/8
    return mpf(alpha*lambda_eff)

###############################
# epsilon factor in Hulthen potential
############################### 

def Hulthen_epsilon(M, z, Isp, nS):
    return kappa * MWTSU2(M, z)/( alpha_eff(Isp, nS) * M )

class Quintuplet_BS:
    """ Set of functions for the Bound States of Quintuplet DM """

    def __init__(self, DarkMatter_Mass, gI, nE, l, Isp, nS, BoundStateFormation_xsec, group_factor):
        self.M = DarkMatter_Mass
        self.gI = gI
        self.nE = nE
        self.l = l
        self.Isp = Isp
        self.nS = nS
        self.bsf = BoundStateFormation_xsec
        self.gf = group_factor

    ###############################
    # Effective couplings of BS
    ###############################

    def alpha_eff(self):
        return alpha_eff(self.Isp, self.nS)

    def h_eps(self, z):
        return Hulthen_epsilon(self.M, z, self.Isp, self.nS)

    # h_eps = Hulthen_epsilon(self.M, z, self.Isp, self.nS)

    # def alpha_eff(self):
    #     lambda_eff = (2 * self.nS**2 - 1 - self.Isp**2)/8
    #     return mpf(alpha*lambda_eff)


    # ###############################
    # # epsilon factor in Hulthen potential
    # ############################### 

    # def Hulthen_epsilon(self, z):
    #     return kappa * MWTSU2(self.M, z)/( self.alpha_eff() * self.M )

    ###############################
    # Bindin energy of BS (for Coulomb potential)
    ###############################

    def binding_energy_BS(self, z):
        prefactor = self.alpha_eff()**2/(4 * self.nE**2)
        H_epsilon = self.h_eps(z)
        corrections = 1 - self.nE**2 * H_epsilon - 0.53 * self.nE**2 * H_epsilon**2 * self.l * (self.l+1)
        return float(prefactor * corrections**2)

    # ###############################
    # # Bindin energy of BS (for Coulomb potential)
    # ###############################

    # def binding_energy_BS(self, z):
    #     prefactor = self.alpha_eff()**2/(4 * self.nE**2)
    #     H_epsilon = self.Hulthen_epsilon(z)
    #     corrections = 1 - self.nE**2 * H_epsilon - 0.53 * self.nE**2 * H_epsilon**2 * self.l * (self.l+1)
    #     return float(prefactor * corrections**2)


    # ###############################
    # # BS number density
    # ###############################


    # def n_BS_eq(self, z):
    #     prefactor = self.gI * (2 * self.M**2/(2 * z * pi))**(3/2)
    #     exponential = np.exp( -(2 - self.binding_energy_BS(z)) * z )
    #     return mpf(prefactor * exponential)


    # ###############################
    # # BS abundance
    # ###############################


    # def Y_BS_eq(self, z):
    #     return self.n_BS_eq(z)/entropy(self.M, z)

    # ###############################
    # # Breaking rates of BS
    # ###############################

    # def Gamma_break(self, z):
    #     pref_1 = gx**2/(2 * self.gI) * self.M**3 * sigma0_prime(self.M)
    #     pref_2 = (1/(z * 4 * pi))**(3/2)
    #     exponential = np.exp( - self.binding_energy_BS(z) * z )
    #     return pref_1 * pref_2 * exponential * self.bsf(z)

    # ###############################
    # # Annihilation rates of BS
    # ###############################

    # def gamma_ann(self):
    #     return self.gf * A25 * self.M

    # def Gamma_ann_Hulthen_TA(self, z):
    #     correction = 1 + self.Hulthen_epsilon(z)**2
    #     return self.gamma_ann() * correction * Kratio(z)

    # def gI_ann_Hulthen(self, z):
    #     prefactor = (2 * self.M**2/(2 * pi * z) )**(3/2)
    #     expon = np.exp( -(2 - a22/2) * z )
    #     return prefactor * expon * self.Gamma_ann_Hulthen_TA(z) 



BS_1s1 = Quintuplet_BS(10 * TeV, 1, 1, 0, 1, 5, BSF_1s1_TA, 3240.0)

# print( BS_1s1.alpha_eff() )
print( BS_1s1.binding_energy_BS(1e6) )
print( BS_1s1.alpha_eff())
# print( BS_1s1.gamma_ann() )

# Gann = [3240.0*A25*mdm, 15625.0/48.0*A25*mdm, 567.0/4.0*A25*mdm, 405.0*A25*mdm, 15625.0/384.0*A25*mdm, 567.0/32.0*A25*mdm, A25*a22*mdm, A25*a22*mdm, A25*a22*mdm ]
# Gdec = [0.0, 0.0, 0.0, s2tw*A25*s2tw*A2MZ*mdm, s2tw*A25*s2tw*A2MZ*mdm, s2tw*A25*s2tw*A2MZ*mdm, 2.0*s2tw*A25*mdm, 1.3*s2tw*A25*mdm, 0.2*s2tw*A25*mdm]
