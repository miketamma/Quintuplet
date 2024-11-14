import numpy as np
from mpmath import *
from scipy import integrate, interpolate
import time
mp.dps = 25; mp.pretty = True
######################################


############################################################################

# Define numerical constants

############################################################################ 

TeV = 1e3

gx = 10 # Quintuplet dof
gstarSS = 106.75 # SM dof
Sqrtgs = sqrt(106.75) # sqrt of SM dof
Mpl = 1.22e19 # Plank mass in GeV
alpha = 0.0313 # alpha EW
a22 = alpha*alpha # alpha^2 EW
A2MZ = 0.0313 # alpha EW
A25 = A2MZ**5.0 # alpha^5 EW (useful for ann gamma)
Tc = 155 # Critical temperature for SM
MW = 80.38 # W mass in GeV
s2tw = 1- (80.38/91.18)**2.0 # sin^2 of weak mixing angle


############################################################################

# SM degrees of freedom as function of temperature

############################################################################ 

Tt, sqrtg = np.loadtxt('Num_lists/sqrtgsT.csv', usecols = (0,1), delimiter = ',', unpack=True)
#Tts, gs = np.loadtxt('dof/gstarsT.csv',usecols = (0,1),delimiter=',',unpack=True)

sqrtg_interp = interpolate.interp1d(Tt, sqrtg, kind='linear')
#interpgs = interp1d(Tts, gs, kind='linear')


############################################################################

# Define entropy and Hubble functions

############################################################################

def entropy(M, z):
    return 2*pi*pi/45 * ( sqrtg_interp(M/z)**2 ) * (M/z)**3

def Hubble(M, z):
    return 2 * np.real( np.sqrt( pi**3/45 * ( sqrtg_interp(M/z)**2 ) ) ) * M**2/(Mpl * z**2)




