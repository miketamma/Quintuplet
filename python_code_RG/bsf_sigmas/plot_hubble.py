import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import matplotlib.tri as tri
from matplotlib.ticker import LogFormatter
from matplotlib.ticker import LogFormatterExponent
import matplotlib.patches as patches
from mpmath import *
from scipy.interpolate import griddata, interpolate, interp1d
import pylab as plab
from labellines import labelLine, labelLines
mpl.style.use('classic')
#####################################

ee,zzi, val1s1,val1s3,val1s5, val2s1,val2s3,val2s5, val2p1,val2p3,val2p5 =np.loadtxt('thavg_bsfall_eps_z.dat',usecols = (0,1,2,3,4,5,6,7,8,9,10),delimiter=',',unpack=True)

Tt, sqrtg  =np.loadtxt('sqrtgsT.csv',usecols = (0,1),delimiter=',',unpack=True)
Tts, gs  =np.loadtxt('gstarsT.csv',usecols = (0,1),delimiter=',',unpack=True)

interpsqgs = interp1d(Tt, sqrtg, kind='linear')
interpgs = interp1d(Tts, gs, kind='linear')

#########################
A2MZ = 0.0313
MV= 80.38
gx2=10.0*10.0
Mpl = 1.22e19
A25 = A2MZ**5.0
s2tw = 1- (80.38/91.18)**2.0

######## Arrays of quanta = 9
# order: 1s1, 1s3, 1s5, 2s1, 2s3, 2s5, 2p1, 2p3, 2p5

lam = [6.0, 5.0, 3.0, 6.0, 5.0, 3.0, 6.0, 5.0, 3.0]
nar = [1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
lar = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
gIr = [1.0, 9.0, 5.0, 1.0, 9.0, 5.0, 3.0, 3.0, 15.0]
###

def ebi(li,nn,ll,mx):
    kappa = 1.74
    yy = kappa * MV/(li*A2MZ*mx)
    pref = (li*A2MZ)**2.0/(4.0*nn*nn)
    term = np.power(1.0 - nn*nn*yy - 0.53*nn*nn*yy*yy*ll*(ll+1.0)  ,2.0)
    return pref*term

def Gbreak(gI,li,nn,ll,mx,zz,sigv):
    s0p = 207.0/20.0*A2MZ**2.0/mx/mx
    ebii = ebi(li,nn,ll,mx)
    term = s0p*sigv*exp(-zz*ebii)
    term1 = (4.0*pi*zz)**(-3.0/2.0)
    pref= gx2/(2.0*gI)*mx*mx*mx
    return pref*term1*term

def lambdadz(zz,mx):
    pref= (45.0/4.0/pi**3.0)**0.5
    term = Mpl/mx/mx
    dof = interpsqgs(mx/zz)/interpgs(mx/zz)
    return pref*term*dof*zz


def lambdabyz2(zz,mx):
    pref= (pi/45.0)**0.5
    term = Mpl*mx
    dof = interpsqgs(mx/zz)
    return pref*term*dof/zz/zz

## build arrays
Gamma1s1=[]
Gamma1s3=[]
Gamma1s5=[]
Gamma2s1=[]
Gamma2s3=[]
Gamma2s5=[]
Gamma2p1=[]
Gamma2p3=[]
Gamma2p5=[]


ldzi=[]
lz2i=[]
mx =13.e3
#mx = np.logspace(3,5,10)
for xx in range(0,len(zzi)):
    # Gamma1s1.append(Gbreak(gIr[0] ,lam[0] ,nar[0] ,lar[0] , mx, zzi[xx], val1s1[xx]))
    # Gamma1s3.append(Gbreak(gIr[1] ,lam[1] ,nar[1] ,lar[1] , mx, zzi[xx], val1s3[xx]))
    # Gamma1s5.append(Gbreak(gIr[2] ,lam[2] ,nar[2] ,lar[2] , mx, zzi[xx], val1s5[xx]))
    # Gamma2s1.append(Gbreak(gIr[3] ,lam[3] ,nar[3] ,lar[3] , mx, zzi[xx], val2s1[xx]))
    # Gamma2s3.append(Gbreak(gIr[4] ,lam[4] ,nar[4] ,lar[4] , mx, zzi[xx], val2s3[xx]))
    # Gamma2s5.append(Gbreak(gIr[5] ,lam[5] ,nar[5] ,lar[5] , mx, zzi[xx], val2s5[xx]))
    # Gamma2p1.append(Gbreak(gIr[6] ,lam[6] ,nar[6] ,lar[6] , mx, zzi[xx], val2p1[xx]))
    # Gamma2p3.append(Gbreak(gIr[7] ,lam[7] ,nar[7] ,lar[7] , mx, zzi[xx], val2p3[xx]))
    # Gamma2p5.append(Gbreak(gIr[8] ,lam[8] ,nar[8] ,lar[8] , mx, zzi[xx], val2p5[xx]))
    ldzi.append(1.0/lambdadz(zzi[xx],mx))
    lz2i.append(1.0/lambdabyz2(zzi[xx],mx))


Gann = [3240.0*A25*mx, 15625.0/48.0*A25*mx, 567.0/4.0*A25*mx, 405.0*A25*mx, 15625.0/384.0*A25*mx, 567.0/32.0*A25*mx, 0.0, 0.0, 0.0 ]
Gdec = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0*s2tw*A25*mx, 1.3*s2tw*A25*mx, 0.2*s2tw*A25*mx]
######## Plot controls
# Define the locations for the axes
left, width = 0.15, 0.80
bottom, height = 0.15, 0.80
bottom_h = left_h = left+width+0.02
    
# Set up the geometry of the three plots
rect_oa = [left, bottom, width, height] # dimensions of temp plot

# Set up the size of the figure
fig = plt.figure(1, figsize=(9.5,9))
# Make the canvas
axoa = plt.axes(rect_oa) 

axoa.set_yscale('log')
axoa.set_xscale('log')
axoa.set_ylim([1.E-24,1.e-9])
axoa.set_xlim([1.e0,1.e6])
axoa.grid(True,alpha=0.5)


axoa.plot(zzi, lz2i, lw=2, ls='-',color='k', label=r'$(\lambda/z^2)^{-1}$')

axoa.plot(zzi, ldzi, lw=2, ls='-.',color='k', label=r'$(\lambda_d \,z)^{-1}$')

axoa.text(1.e2,1.e-21,r'Hubble coefficients', fontsize=20)
axoa.text(1.e2,1.e-22,r' $m_\chi =$'+str(mx/1.e3)+' TeV', fontsize=20)


#handles,labels = axoa.get_legend_handles_labels()

axoa.minorticks_on()
axoa.tick_params('both', length=5, width=2, which='major',labelsize=15)
axoa.tick_params('both', length=4.5, width=2, which='minor',labelsize=13)
axoa.tick_params('x', length=5, width=2, which='major',labelsize=18)


axoa.set_xlabel(r'$z$ ',fontsize=25)

axoa.set_ylabel(r'$(\lambda/z^2)^{-1}\,[{\rm GeV^{-2}}]\quad {\rm and}\quad  (\lambda_d z)^{-1}\,[{\rm GeV}]$',fontsize=25)


p=patches.Rectangle((0, 0), 1, 1,fill=False, transform=axoa.transAxes, clip_on=False,lw=2)
axoa.add_patch(p) 

#leg=axoa.legend(handles,labels,bbox_to_anchor=(0.25, 0.95), loc=2, borderaxespad=0.,fancybox=True, framealpha=0.1)
#leg.get_frame().set_linewidth(0.0)
labelLines(axoa.get_lines(), xvals=(1, 100), zorder=2.5,fontsize =25)

plab.savefig('hubble.pdf', bbox_inches=0,dpi=100)

plt.show()
