import numpy as np
from mpmath import *
from scipy import integrate, interpolate
import time
mp.dps = 25; mp.pretty = True
######################################

start = time.time()

g_chi = 10

#### Define alpha_2 and alpha = pi*alpha_2 #####

A2MZ = 0.0313
# alpha = mpf(pi*A2MZ)

def alpha_eff(Isp, n):
    lambda_eff = (2*n*n - 1 - Isp*Isp)/8
    return mpf(pi*A2MZ*lambda_eff)


#### Import high precision Bessel function tables and define the interpolation #####

K1_hiprec = np.loadtxt('K1_hiprec.txt')
K2_hiprec = np.loadtxt('K2_hiprec.txt')

k1_interp = interpolate.interp1d( K1_hiprec[:, 0], K1_hiprec[:, 1] )
k2_interp = interpolate.interp1d( K2_hiprec[:, 0], K2_hiprec[:, 1] )

def K1(z):
    return 10**mpf( float( k1_interp(z) ) )

def K2(z):
    return 10**mpf( float( k2_interp(z) ) )


##### Define Sommerfeld effect with Huelten potential ####

def SommerfeldHuelten(XX, EE):
    kaps= 1.74
    if(EE==0):
        res= 2.0*XX/(1-exp(-2.0*XX))
    else:
        arg1 = (pi*pi)*kaps/fabs(XX*EE)
        arg2 = sqrt( 1.0 - ( 4.0/(pi*pi*kaps) )*XX*XX*EE ) 
        num = 2*XX*sinh( arg1 )
        den = cosh( arg1 ) - cosh( arg1*arg2 )
        res= num/den
    return res

def NoSommerfeld(X, E):
    return 1

##### Define Thermal Average integrand with any function of (X, eps) ####

def ThermalAverage(Sfunc, Isp, n, beta, eps, z):

    X = alpha_eff(Isp,n)/beta

    den = 2*K2(z)*K2(z)
    gamma = float( 1.0/sqrt(1.0 - beta*beta/4.0) )
    t2 =  (beta**2)*(gamma**7)*K1(2.0*z*gamma)
    num = Sfunc(X, eps)*t2
    res = z*num/den
    return res

## Define a grid in epsilon and z and do the thermal average #####

nQuint = 5 # Fix Quintuplet
limits = np.logspace(-9, np.log10(1.999), 100)

zarr = np.logspace(0, 6.0, 500)
EEarr = np.logspace(-5, 2.0, 300)

### Do Thermal Average for NoSommerfeld ###

ListTA = []

for zi in range(0, len(zarr)):
    for eei in range(0, len(EEarr)):

        thav_S1_list = []
        thav_S3_list = []
        thav_S5_list = [] 
        for beta in limits:
            thav_S1_list.append( ThermalAverage(NoSommerfeld, 1, nQuint, beta, EEarr[eei], zarr[zi]) )
            thav_S3_list.append( ThermalAverage(NoSommerfeld, 3, nQuint, beta, EEarr[eei], zarr[zi]) )
            thav_S5_list.append( ThermalAverage(NoSommerfeld, 5, nQuint, beta, EEarr[eei], zarr[zi]) )
        val_S1 = mpf( re( integrate.trapezoid(thav_S1_list, limits) ) )
        val_S3 = mpf( re( integrate.trapezoid(thav_S3_list, limits) ) )
        val_S5 = mpf( re( integrate.trapezoid(thav_S5_list, limits) ) )
        ListTA.append(np.asarray( [EEarr[eei], zarr[zi], sign(val_S1)*mpf(log10(fabs(val_S1))), sign(val_S3)*mpf(log10(fabs(val_S3))), sign(val_S5)*mpf(log10(fabs(val_S5)))]) )

np.savetxt('NoSommTA.txt', ListTA)


## Do Themarl Average for SommerfeldHuelten ###

ListTA = []

for zi in range(0, len(zarr)):
    for eei in range(0, len(EEarr)):

        limits = np.logspace(-9, np.log10(0.999), 100)
        thav_S1_list = []
        thav_S3_list = []
        thav_S5_list = [] 
        for beta in limits:
            thav_S1_list.append( ThermalAverage(SommerfeldHuelten, 1, nQuint, beta, EEarr[eei], zarr[zi]) )
            thav_S3_list.append( ThermalAverage(SommerfeldHuelten, 3, nQuint, beta, EEarr[eei], zarr[zi]) )
            thav_S5_list.append( ThermalAverage(SommerfeldHuelten, 5, nQuint, beta, EEarr[eei], zarr[zi]) )
        val_S1 = mpf( re( integrate.trapezoid(thav_S1_list, limits) ) )
        val_S3 = mpf( re( integrate.trapezoid(thav_S3_list, limits) ) )
        val_S5 = mpf( re( integrate.trapezoid(thav_S5_list, limits) ) )
        ListTA.append(np.asarray( [EEarr[eei], zarr[zi], sign(val_S1)*mpf(log10(fabs(val_S1))), sign(val_S3)*mpf(log10(fabs(val_S3))), sign(val_S5)*mpf(log10(fabs(val_S5)))]) )


np.savetxt('SommHueltenTA.txt', ListTA)

end = time.time()
print(end - start)





