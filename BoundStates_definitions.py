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


# ###############################
# # Effective couplings of BS
# ###############################

# def alpha_eff(Isp, nS):
#     lambda_eff = (2 * nS**2 - 1 - Isp**2)/8
#     return mpf(alpha*lambda_eff)

# ###############################
# # epsilon factor in Hulthen potential
# ############################### 

# def Hulthen_epsilon(M, z, Isp, nS):
#     return kappa * MWTSU2(M, z)/( alpha_eff(Isp, nS) * M )

class Quintuplet_BS(Quintuplet_DM):
    """Properties and reaction rates of a quintuplet-DM bound state."""

    def __init__(self, DarkMatter_Mass, gI, nE, l, Isp, nS, BoundStateFormation_xsec, group_factor):
        """Store the bound-state quantum numbers and tabulated BSF rate."""
        super().__init__(DarkMatter_Mass)
        self.gI = gI
        self.nE = nE
        self.l = l
        self.Isp = Isp
        self.nS = nS
        self._bsf_table = BoundStateFormation_xsec
        self.gf = group_factor
        Quintuplet_DM.__init__(self, DarkMatter_Mass)

    ###############################
    # Effective couplings of BS
    ###############################

    def lambda_eff(self):
        """Return the attractive-potential group coefficient."""
        return (2 * self.nS**2 - 1 - self.Isp**2) / 8

    def alpha_eff_BSF(self):
        """Return the effective coupling at the bound-state formation scale."""
        scale = self.lambda_eff() * alpha2(MZ) * self.M
        return self.lambda_eff() * alpha2(scale)

    def alpha_eff_BE(self):
        """Return the effective coupling at the binding-energy scale."""
        scale = (self.lambda_eff() * alpha2(MZ))**2 * self.M
        return self.lambda_eff() * alpha2(scale)


    ###############################
    # Screening in Hulthen potential
    ############################### 

    def Hulthen_epsilon(self, alpha, z):
        """Return the dimensionless screening parameter of the Hulthén potential."""
        return kappa * MWTSU2(self.M, z)/( alpha() * self.M )

    def screening_correction(self, z):
        """Return the finite-range correction entering the binding energy."""
        epsilon = self.Hulthen_epsilon(self.alpha_eff_BE, z)

        return 1 - self.nE**2 * epsilon - 0.53 * self.nE**2 * epsilon**2 * self.l * (self.l+1)

    def is_bound(self, z):
        """Return whether the screened potential supports this bound state."""
        return self.screening_correction(z) > 0


    ###############################
    # Binding energy of BS (for Coulomb potential)
    ###############################

    def binding_energy_BS(self, z):
        """Return the binding energy divided by the constituent DM mass."""
        correction = self.screening_correction(z)

        # The state has dissolved when the screening correction is non-positive.
        if correction <= 0:
            return 0.0

        coulomb_energy = self.alpha_eff_BE()**2 / (4 * self.nE**2)
        return float(coulomb_energy * correction**2)

    def bsf(self, z):
        """Return the tabulated BSF factor, or zero if the state is unbound."""
        if not self.is_bound(z):
            return 0.0

        return self._bsf_table(z)



    def n_BS_eq(self, z):
        """Return the equilibrium number density of the bound state."""
        if not self.is_bound(z):
            return mpf("0.0")
        prefactor = self.gI * (2 * self.M**2/(2 * z * pi))**(3/2)
        exponential = np.exp( -(2 - self.binding_energy_BS(z)) * z )
        return mpf(prefactor * exponential)

    def Y_BS_eq(self, z):
        """Return the equilibrium bound-state yield n_BS/s."""
        return self.n_BS_eq(z) / self.entropy(z)

    def Gamma_break(self, z):
        """Return the thermally suppressed bound-state dissociation rate."""
        if not self.is_bound(z):
            return 0.0

        scale = self.lambda_eff() * alpha2(MZ) * self.M
        prefactor = gx**2/(2 * self.gI) * self.M**3 * sigma0_prime( scale, self.M)
        thermal_factor = (1/(z * 4 * pi))**(3/2)
        exponential = np.exp( - self.binding_energy_BS(z) * z )
        return prefactor * thermal_factor * exponential * self.bsf(z) 

    def Gamma_break_NoExp(self, z):
        """Return the dissociation rate without its Boltzmann exponential."""
        if not self.is_bound(z):
            return 0.0

        scale = self.lambda_eff() * alpha2(MZ) * self.M

        prefactor = gx**2/(2 * self.gI) * self.M**3 * sigma0_prime(( self.lambda_eff(self) * alpha2(MZ) ) * self.M, self.M)
        thermal_factor = (1/(z * 4 * pi))**(3/2)
        return prefactor * thermal_factor * self.bsf(z)

    ###############################
    # Annihilation rates of BS
    ###############################

    def gamma_ann(self):
        """Return the dimensionless short-distance annihilation coefficient."""
        return self.gf * A25

    def Gamma_ann_Hulthen_TA(self, z):
        """Return the thermally averaged bound-state annihilation rate."""
        if not self.is_bound(z):
            return 0.0

        epsilon = self.Hulthen_epsilon( self.alpha_eff_BE, z)
        hulthen_correction = 1 + epsilon**2
        
        return self.gamma_ann() * self.M * hulthen_correction * Kratio(z)

    def gI_ann_Hulthen(self, z):
        """Return the equilibrium annihilation reaction density."""
        if not self.is_bound(z):
            return 0.0
        prefactor = (2 * self.M**2/(2 * pi * z) )**(3/2)
        exponential = np.exp( -(2 - a22/2) * z )
        return prefactor * exponential * self.Gamma_ann_Hulthen_TA(z) 




# Variable order is: (Mass, gI, nE, l, I, nS, BS function, group factor constant)

# BS order is: (1s1, 1s3, 1s5, 2s1, 2s3, 2s5, 2p1, 2p3, 2p5)

nS_Quint = 5.0

gI_list = [1.0, 9.0, 5.0, 1.0, 9.0, 5.0, 3.0, 3.0, 15.0]
nE_list = [1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
l_list = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
Isp_list = [1.0, 3.0, 5.0, 1.0, 3.0, 5.0, 1.0, 3.0, 5.0]

BS_func_list = [BSF_1s1_TA, BSF_1s3_TA, BSF_1s5_TA, BSF_2s1_TA, BSF_2s3_TA, BSF_2s5_TA, BSF_2p1_TA, BSF_2p3_TA, BSF_2p5_TA]

gf_list = [3240.0, 15625.0/48.0, 567.0/4.0, 405.0, 15625.0/384.0, 567.0/32.0, 1, 1, 1]

# Gdec = [0.0, 0.0, 0.0, s2tw*A25*s2tw*A2MZ*mdm, s2tw*A25*s2tw*A2MZ*mdm, s2tw*A25*s2tw*A2MZ*mdm, 2.0*s2tw*A25*mdm, 1.3*s2tw*A25*mdm, 0.2*s2tw*A25*mdm]

# BS_1s1 = Quintuplet_BS(10 * TeV, gI_list[0], nE_list[0], l_list[0], Isp_list[0], nS_Quint, BS_func_list[0], gf_list[0])

# print(BS_1s1.gI)

# print( BS_1s1.entropy(1e3) )
# print( BS_1s1.Y_BS_eq(1e3) )
