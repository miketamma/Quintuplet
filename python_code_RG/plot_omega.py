import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import matplotlib.tri as tri
from matplotlib.ticker import LogFormatter
from matplotlib.ticker import LogFormatterExponent
import matplotlib.patches as patches
import pylab as plab
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

mpl.style.use('classic')
#####################################

mdm, omega_3bseff, omega_3bseffgsc, omega_net, omega_netibs=np.loadtxt('mx_omega_huelthen_network.dat',usecols = (0,1,2,3,4),delimiter=',',unpack=True)



#######################################
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


#axoa.set_yscale('log')
#axoa.set_xscale('log')
axoa.set_ylim([1.E-3,0.2])
axoa.set_xlim([2.e0,15])
axoa.grid(True,alpha=0.5)


inset_ax = inset_axes(axoa, width="50%", height=2.5,loc=2)
inset_ax.set_xlim([12,14])
inset_ax.set_ylim([0.11,0.13])
inset_ax.grid(True,alpha=0.5)


fac = 1.e3

axoa.plot(mdm/fac, omega_3bseffgsc, c='g', label=r'Effective: const. $g_\star$',lw=3 , ls=':')

axoa.plot(mdm/fac, omega_3bseff, c='g', label=r'Effective: $g_\star (T)$',lw=3 , ls='-.')

axoa.plot(mdm/fac, omega_net, c='k', label=r'Network: $g_\star (T)$',lw=3 , ls='--')

axoa.plot(mdm/fac, omega_netibs, c='k', label=r'Network: $g_\star (T)$ + $B_I \rightarrow B_J$',lw=3 , ls='-')


inset_ax.plot(mdm/fac, omega_3bseffgsc, c='g', label=r'Effective: const. $g_\star$',lw=3 , ls=':')

inset_ax.plot(mdm/fac, omega_3bseff, c='g', label=r'Effective: $g_\star (T)$',lw=3 , ls='-.')

inset_ax.plot(mdm/fac, omega_net, c='k', label=r'Network: $g_\star (T)$',lw=3 , ls='--')

inset_ax.plot(mdm/fac, omega_netibs, c='k', label=r'Network: $g_\star (T)$ + $B_I \rightarrow B_J$',lw=3 , ls='-')

inset_ax.axhline(y=0.119,color='r',lw = 2)


axoa.axhline(y=0.119,color='k',lw = 2)
axoa.text(7.5,0.02,'Huelthen Potential + massless B-S',fontsize=15)


axoa.set_xlabel(r'$m_\chi\,[{\rm TeV}]$ ',fontsize=25)

axoa.set_ylabel(r'$\Omega h^2$',fontsize=25)


handles,labels = axoa.get_legend_handles_labels()

axoa.minorticks_on()
axoa.tick_params('both', length=5, width=2, which='major',labelsize=15)
axoa.tick_params('both', length=4.5, width=2, which='minor',labelsize=13)
axoa.tick_params('x', length=5, width=2, which='major',labelsize=18)


p=patches.Rectangle((0, 0), 1, 1,fill=False, transform=axoa.transAxes, clip_on=False,lw=2)
axoa.add_patch(p) 

leg=axoa.legend(handles,labels,bbox_to_anchor=(0.05, 0.45), loc=2, borderaxespad=0.,fancybox=True, framealpha=0.1)
leg.get_frame().set_linewidth(0.0)

plab.savefig('omega_compare_all.pdf', bbox_inches=0,dpi=100)

plt.show()


