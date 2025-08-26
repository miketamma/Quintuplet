import numpy as np
from mpmath import *
from scipy import integrate, interpolate
from scipy.integrate import solve_ivp
import time
mp.dps = 25; mp.pretty = True
######################################


############################################################################

# Define numerical constants

############################################################################ 

eV = 1e-9
TeV = 1e3

gx = 10 # Quintuplet dof
gstarSS = 106.75 # SM dof
Sqrtgs = sqrt(106.75) # sqrt of SM dof
Mpl = 1.22e19 # Plank mass in GeV
alpha = 0.0313 # alpha EW
a22 = alpha*alpha # alpha^2 EW
A2MZ = 0.0313 # alpha EW
A25 = A2MZ**5.0 # alpha^5 EW (useful for ann gamma)
Tc = 155.0 # Critical temperature for SM
MW = 80.38 # W mass in GeV
s2tw = 1- (80.38/91.18)**2.0 # sin^2 of weak mixing angle

kappa = 1.74 # Hulthen factor 

nS_Quintuplet = 5 # Quintuplet state


############################################################################

# SM degrees of freedom as function of temperature

############################################################################ 

Tt, sqrtg = np.loadtxt('Num_lists/sqrtgsT.csv', usecols = (0,1), delimiter = ',', unpack=True)
Tts, gs = np.loadtxt('Num_lists/gstarsT.csv',usecols = (0,1),delimiter=',',unpack=True)

sqrt_gs_eff = interpolate.interp1d(Tt, sqrtg, kind='linear')
gs_star = interpolate.interp1d(Tts, gs, kind='linear')


def sqrt_gs(T):
    return gs_star(T) / sqrt_gs_eff(T)

############################################################################

# Define entropy, Hubble and DM energy density functions

############################################################################


def Omega_DM(h, Y_DM_inf, M):
    return (0.11/h**2) * (Y_DM_inf * M)/(0.4 * eV)


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
# Definie ratio K1/K2
###############################

def Kratio(z):
    return K1(z)/K2(z)



############################################################################

# Definie numerical Boltzmann ODE solver, n is the number of BS variables

############################################################################

def NDE_solver(Boltz, n):
    y0 = np.concatenate( ([1.e-30], np.full( n, 1.e-30) ) ) #np.full( n, 1.e-30)
    teval = np.logspace(np.log10(3.5), 6.9, 100)
    tspan=[3.0, 1.e7]
    
    return solve_ivp(Boltz, tspan, y0 = y0, method = 'Radau', t_eval = teval, atol=1.e-20)