# Computation of the Quintuplet relic thermal mass

Precise calculation of the relic mass of the Dark Matter, assumed as being part of the 5 representation of the Standard Model SU(2) gauge group.


## TL;DR

The final result can be obtained by simply running the 'Boltzmann_solver.py' script; it will ask as input the number of mass points to scan in the M$_{DM}$ ∈ [3,15] TeV interval, and which case to run (more details below). It outputs the list [M$_{DM}$, Ω$_{DM}$ h$^{2}$] as a .txt file and optionally it prints it on terminal.

## Structure of the code

The code is structured in a subset of scripts that contain the different definitions necessary to compute the Boltzmann equations.


### Num_lists

This folder contains all the necessary pretabulated quantities as .txt files. Scripts to (re)compute these tables are provided.

- 'K1_hiprec.txt' and 'K2_hiprec.txt' contain the modified Bessel functions of first and second kind, K$_{1}$ and K$_{2}$, as function of 
-
-

In Num_lists folder there are the txt files for K1 and K2 (Bessel functions), computed with high precision, and the thermally averaged cross sections for s-wave production with Hulthen potential and bound state formation (BSF) for the states 1s1, 1s32 1s5, 2s1, 2s3, 2s5, 2p1, 2p3, 2p5

Basic_definitions.py - contains constant definitions and import Bessel functions

Quintuplet_definitions.py - contains the Quintuplet_DM class definition

BoundStates_definitions.py - contains the Quintuplet_BS class definition

NDE_solver.py - contains the ODE solver function, and does the mass scan for different cases (IMPLEMENTED: Hulthen_free, Hulthen_1s_effective) 
