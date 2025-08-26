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


########### Define ArcCotan with numpy

def arcctg(x):
    return np.arctan(1/x)

##### Define bound state formation cross sections ####


##### 1s (n = 1, l = 0) ####

def sv_10(s, li, lf, CJ, Ctau, XX):

    pref = li*( (lf*XX)**5 ) * (2*s + 1)/( g_chi**2 )
    num = ( 2**11 ) * pi * ( 1 + XX*XX*li*li ) * exp(-4*XX*li*arcctg(XX*lf))
    den = 3*( (1 + XX*XX*lf*lf)**3 )*(1 - exp(-2*pi*XX*li) )
    group = (CJ + Ctau/lf)**2

    return pref * num * group/den

##### 1s1 ####

def sv_1s1_bsf(XX):

    return sv_10(0, 5, 6, sqrt(6), -sqrt(6), XX)

##### 1s3 ####

def sv_1s3_bsf(XX):
    return sv_10(1, 6, 5, sqrt(6), sqrt(6), XX) + sv_10(1, 3, 5, sqrt(21)/sqrt(2), - sqrt(42), XX)


##### 1s5 ####

def sv_1s5_5_bsf(XX):

    CJ, Ctau = sqrt(12), -3*sqrt(12)
    
    num = 2304*( (3*CJ + Ctau)**2 )*(XX**4)
    den = 25*( (1 + 9*XX*XX)**3 )

    return num/den 

def sv_1s5_bsf(XX):
    return sv_10(0, 5, 3, sqrt(21)/sqrt(2), sqrt(42), XX) + sv_1s5_5_bsf(XX)


##### 2s (n = 2, l = 0) ####

def sv_20(s, li, lf, CJ, Ctau, XX):

    pref = li * ( (lf*XX)**5 ) * (2*s + 1)/( g_chi**2 )
    num = ( 2**14 ) * pi * ( 1 + XX*XX*li*li ) * exp(-4*XX*li*arcctg(XX*lf/2))
    den = 3 * ( (4 + XX*XX*lf*lf)**5 ) * (1 - exp(-2*pi*XX*li) )
    group = ( CJ * ( XX*XX*lf * ( lf - 2 * li) - 4 ) + Ctau * ( XX*XX*(3*lf - 4*li) - 4/lf ) )**2

    return pref * num * group/den

##### 2s1 ####

def sv_2s1_bsf(XX):

    return sv_20(0, 5, 6, sqrt(6), -sqrt(6), XX)


##### 2s3 ####

def sv_2s3_bsf(XX):
    return sv_20(1, 6, 5, sqrt(6), sqrt(6), XX) + sv_20(1, 3, 5, sqrt(21)/sqrt(2), - sqrt(42), XX)


##### 2s5 ####

def sv_2s5_5_bsf(XX):

    CJ, Ctau = sqrt(12), -3*sqrt(12)
    
    num = 18432*( (-12*CJ  -4*Ctau + 27*(CJ + Ctau)*XX*XX)**2 )*(XX**4)
    den = 25*( (4 + 9*XX*XX)**5 ) 

    return num/den

def sv_2s5_bsf(XX):
    return sv_20(0, 5, 3, sqrt(21)/sqrt(2), sqrt(42), XX) + sv_2s5_5_bsf(XX)


##### 2p (n = 2, l = 1) ####

def sv_21(s, li, lf, CJ, Ctau, XX):

    pref = li*( (lf*XX)**5 ) * (2*s + 1)/( g_chi**2 )
    num = ( 2**12 ) * pi * A2MZ * ( XX**2 ) * exp(-4*XX*li*arcctg(XX*lf/2))
    den = 9 * ( (4 + XX*XX*lf*lf)**5 )*(1 - exp(-2*pi*XX*li) )
    group1 = ( CJ * ( lf * ( XX*XX*li*( 3*lf - 4*li) + 8 ) - 12 * li ) + Ctau * ( XX*XX * ( -3*lf*lf + 12*li*lf - 8*li*li) + 4 ) )**2
    group2 = ( 2**5 ) * ( XX*XX*li*li + 1 ) * ( XX*XX*li*li + 4 ) * ( CJ * lf + 2 * Ctau )**2

    return pref * num * ( group1 + group2 )/den

##### 2p1 ####

def sv_2p1_bsf(XX):
    return sv_21(1, 5, 6, sqrt(6), -sqrt(6), XX) + sv_21(1, 3, 6, sqrt(6), -sqrt(6), XX)


##### 2p3 ####

def sv_2p3_bsf(XX):
    return 2*sv_21(0, 6, 5, sqrt(6), sqrt(6), XX) + 2 * sv_21(0, 3, 5, sqrt(21)/sqrt(2), -sqrt(42), XX)

##### 2p5 ####

def sv_2p5_2_bsf(XX):

    CJ, Ctau = sqrt(12), -3*sqrt(12)
    
    num = 41472*A2MZ*( 128*( (3*CJ + 2*Ctau)**2 ) + ( 24*CJ + Ctau*(4 - 27*XX*XX) )**2 )*(XX**6)
    den = 25*( (4 + 9*XX*XX)**5 ) 

    return num/den

def sv_2p5_bsf(XX):
    return 2 * sv_21(1, 5, 3, sqrt(21)/sqrt(2), sqrt(42), XX) + 2 * sv_2p5_2_bsf(XX)

##### Define Thermal Average integrand with any function of (X, eps) ####

def ThermalAverage_bsf_Massless(Sfunc, beta, z):

    X = A2MZ/beta

    den = 2*K2(z)*K2(z)
    gamma = float( 1.0/sqrt(1.0 - beta*beta/4.0) )
    t2 =  (beta**2)*(gamma**7)*K1(2.0*z*gamma)
    num = Sfunc(X)*t2
    res = z*num/den
    return res


## Define a grid in epsilon and z and do the thermal average #####
nQuint = 5 # Fix Quintuplet
limits = np.logspace(-9, np.log10(1.999), 100)

zarr = np.logspace(0, 7.0, 700)


# Do Thermal Average for Bound States ###


#### 1s ####

ListTA = []

for zi in range(0, len(zarr)):

    thav_list = []
    for beta in limits:
        thav_list.append( ThermalAverage_bsf_Massless(sv_1s1_bsf, beta, zarr[zi]) )
    val = mpf( re( integrate.trapezoid(thav_list, limits) ) )
    ListTA.append(np.asarray( [zarr[zi], sign(val)*mpf(log10(fabs(val))) ]) )

np.savetxt('./BSF_1s1_TA.txt', ListTA)

ListTA = []

for zi in range(0, len(zarr)):

    thav_list = []
    for beta in limits:
        thav_list.append( ThermalAverage_bsf_Massless(sv_1s3_bsf, beta, zarr[zi]) )
    val = mpf( re( integrate.trapezoid(thav_list, limits) ) )
    ListTA.append(np.asarray( [zarr[zi], sign(val)*mpf(log10(fabs(val))) ]) )

np.savetxt('./BSF_1s3_TA.txt', ListTA)

ListTA = []

for zi in range(0, len(zarr)):

    thav_list = []
    for beta in limits:
        thav_list.append( ThermalAverage_bsf_Massless(sv_1s5_bsf, beta, zarr[zi]) )
    val = mpf( re( integrate.trapezoid(thav_list, limits) ) )
    ListTA.append(np.asarray( [zarr[zi], sign(val)*mpf(log10(fabs(val))) ]) )

np.savetxt('./BSF_1s5_TA.txt', ListTA)


##### 2s ####


ListTA = []

for zi in range(0, len(zarr)):

    thav_list = []
    for beta in limits:
        thav_list.append( ThermalAverage_bsf_Massless(sv_2s1_bsf, beta, zarr[zi]) )
    val = mpf( re( integrate.trapezoid(thav_list, limits) ) )
    ListTA.append(np.asarray( [zarr[zi], sign(val)*mpf(log10(fabs(val))) ]) )

np.savetxt('./BSF_2s1_TA.txt', ListTA)

ListTA = []

for zi in range(0, len(zarr)):

    thav_list = []
    for beta in limits:
        thav_list.append( ThermalAverage_bsf_Massless(sv_2s3_bsf, beta, zarr[zi]) )
    val = mpf( re( integrate.trapezoid(thav_list, limits) ) )
    ListTA.append(np.asarray( [zarr[zi], sign(val)*mpf(log10(fabs(val))) ]) )

np.savetxt('./BSF_2s3_TA.txt', ListTA)

ListTA = []

for zi in range(0, len(zarr)):

    thav_list = []
    for beta in limits:
        thav_list.append( ThermalAverage_bsf_Massless(sv_2s5_bsf, beta, zarr[zi]) )
    val = mpf( re( integrate.trapezoid(thav_list, limits) ) )
    ListTA.append(np.asarray( [zarr[zi], sign(val)*mpf(log10(fabs(val))) ]) )

np.savetxt('./BSF_2s5_TA.txt', ListTA)


##### 2p ####


ListTA = []

for zi in range(0, len(zarr)):

    thav_list = []
    for beta in limits:
        thav_list.append( ThermalAverage_bsf_Massless(sv_2p1_bsf, beta, zarr[zi]) )
    val = mpf( re( integrate.trapezoid(thav_list, limits) ) )
    ListTA.append(np.asarray( [zarr[zi], sign(val)*mpf(log10(fabs(val))) ]) )

np.savetxt('./BSF_2p1_TA.txt', ListTA)

ListTA = []

for zi in range(0, len(zarr)):

    thav_list = []
    for beta in limits:
        thav_list.append( ThermalAverage_bsf_Massless(sv_2p3_bsf, beta, zarr[zi]) )
    val = mpf( re( integrate.trapezoid(thav_list, limits) ) )
    ListTA.append(np.asarray( [zarr[zi], sign(val)*mpf(log10(fabs(val))) ]) )

np.savetxt('./BSF_2p3_TA.txt', ListTA)

ListTA = []

for zi in range(0, len(zarr)):

    thav_list = []
    for beta in limits:
        thav_list.append( ThermalAverage_bsf_Massless(sv_2p5_bsf, beta, zarr[zi]) )
    val = mpf( re( integrate.trapezoid(thav_list, limits) ) )
    ListTA.append(np.asarray( [zarr[zi], sign(val)*mpf(log10(fabs(val))) ]) )

np.savetxt('./BSF_2p5_TA.txt', ListTA)

end = time.time()
print(end - start)





