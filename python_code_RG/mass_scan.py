import matplotlib.pyplot as plt
import numpy as np
import os, sys, subprocess, glob, time, random, shutil
import time
############################
t0 = time.time()  # start timing


massarr=np.logspace(np.log10(2.e3), np.log10(15.e3) ,20)
fileout = open("mx_omega_huelthen_network_3.dat", "w")


for ii in range(0,len(massarr)):
    cmd="python3 boltzmann_eqs_solve.py"+"  "+str(massarr[ii])
    subprocess.call(cmd, shell=True)
    mdm, omega_3bseff, omega_3bseffgsc, omega_net, omega_netibs=np.loadtxt('tempr.dat',usecols = (0,1,2,3,4),delimiter=',',unpack=True)
    fileout.write("%e \t ,  %e \t ,  %e \t ,  %e \t ,  %e   \n" %(mdm, omega_3bseff, omega_3bseffgsc, omega_net, omega_netibs ) )

fileout.close()

tf = time.time()

print('Total time taken == ', tf-t0)
