import matplotlib.pyplot as plt
import numpy as np
from mpmath import *
mp.dps = 30; mp.pretty = True
######################################

def trapz(tabx, taby):
    res=0
    for pp in range(0,len(tabx)):
        if(pp==0):
            res= 0.5*taby[pp]*tabx[pp]
        else:
            res= res + 0.5*(taby[pp]+ taby[pp-1] )*(tabx[pp]-tabx[pp-1])
    return res

def Somef(XX, EE):
    pi2 = pi*pi
    kaps= 1.74
    if(EE==0):
        res= 2.0*XX/(1-exp(-2.0*XX))
    else:
        res= 2*XX*sinh(pi2/fabs(XX*EE*kaps))/(cosh(pi2/fabs(XX*EE*kaps)) - cosh(pi2/(fabs(XX*EE*kaps))*sqrt(1.0+ 4.0/(pi*pi)*XX*XX*EE*kaps)))
    return res


def thavgfunc(BB,AA,EE,ZZ):
    pref = ZZ/2.0
    den =  besselk(2, ZZ)*besselk(2, ZZ)
    XX = AA/BB
    t1 = BB*BB
    t2 = (4/(4-BB*BB))**(7/2)
    t3 = 4/sqrt(4-BB*BB)
    num= Somef(XX, EE)*t1*t2*besselk(1, ZZ*t3)
    res = pref*num/den
    return res

A2MZ = 0.0313

zarr = np.logspace(np.log10(3.0),6.0,200)
#zarr = np.logspace(np.log10(25),6.0,50)
EEarr1 = np.logspace(-3,1.0,1000)
newarr = -1.0*EEarr1
new0= [0]
EEarr = [*newarr,*new0]
fileout = open("thavg_huelthen_eps_z_s653.dat", "w")

count=0

for zi in range(0,len(zarr)):
    for eei in range(0,len(EEarr)):
        #for zi in range(5,6):
        #    for eei in range(11,12):
        argee= mpf(EEarr[eei])
        argzz = mpf(zarr[zi])
        alpah = mpf(pi*A2MZ)
        alpah6 = mpf(6*pi*A2MZ)
        alpah5 = mpf(5*pi*A2MZ)
        alpah3 = mpf(3*pi*A2MZ)
        limits= [mpf(1.e-5),mpf(0.5)]
        #print('inputs eps, z =  ',argee,'  ,  ',argzz)
        barr=np.logspace(-9,np.log10(1.999),100)
        bbtab=[]
        thfunctab=[]
        thfunctab6=[]
        thfunctab5=[]
        thfunctab3=[]
        t6=0
        t5=0
        t3=0
        for kk in range(0, len(barr)):
            bbtab.append(mpf(barr[kk]))
            t6= re(thavgfunc(mpf(barr[kk]),alpah6,argee,argzz))
            t5= re(thavgfunc(mpf(barr[kk]),alpah5,argee,argzz))
            t3= re(thavgfunc(mpf(barr[kk]),alpah3,argee,argzz))
            thfunctab6.append(t6)
            thfunctab5.append(t5)
            thfunctab3.append(t3)
            thfunctab.append(t6+t5+t3)
        val = fabs(trapz(bbtab, thfunctab))
        val6 = fabs(trapz(bbtab, thfunctab6))
        val5 = fabs(trapz(bbtab, thfunctab5))
        val3 = fabs(trapz(bbtab, thfunctab3))
        count=count+1
        fileout.write("%e \t ,  %e \t ,  %e \t ,  %e \t ,  %e \t ,  %e \t ,  %e \t ,  %e \t ,  %e \t ,  %e  \n" %(argee,argzz, val, log10(val), val6, log10(val6), val5, log10(val5), val3, log10(val3)  ) )
        # val,err = quad(lambda x: thavgfunc(x,alpah,argee,argzz),limits,maxdegree= 10,error=True)
        # res= fabs(re(val))
        # count=count+1
        # fileout.write("%e \t  %e \t  %e \t  %e \t  %e \t  %e \n" %(argee,argzz, res,err, log10(res),log10(err)) )
        print(argee,'  ,  ',argzz, '  ,  ', val, '  ,  ', log10(val) , '  ,  ', val6, '  ,  ',log10(val6) )

fileout.close()

