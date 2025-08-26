import matplotlib.pyplot as plt
import numpy as np
from mpmath import *
mp.dps = 25; mp.pretty = True

xs = np.logspace(0, 9.0, 30000)


k1_hiprec = []
k2_hiprec = []

nDMeq_z = []
YDMeq_z = []

# print( mpf( log10( besselk(1, 10**6) ) ) ) 

for x in range(0,len(xs)):
    k1_hiprec.append( np.asarray( [ float(xs[x]), log10( besselk(1, xs[x])) ] ) )
    k2_hiprec.append( np.asarray( [ float(xs[x]), log10( besselk(2, xs[x])) ] ) )
    nDMeq_z.append( np.asarray( [ float(xs[x]), log10( besselk(2, xs[x])/xs[x]) ] ) )
    YDMeq_z.append( np.asarray( [ float(xs[x]), log10( besselk(2, xs[x])*xs[x]*xs[x]) ] ) )


np.savetxt('./K1_hiprec.txt', k1_hiprec)
np.savetxt('./K2_hiprec.txt', k2_hiprec)
np.savetxt('./nDMeq_z.txt', nDMeq_z)
np.savetxt('./YDMeq_z.txt', YDMeq_z)
