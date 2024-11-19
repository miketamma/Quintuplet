# Quintuplet
Precise calculation of the mass of Dark Matter, as Quintuplet under SU(2).

In Num_lists folder there are the txt files for K1 and K2 (Bessel functions), computed with high precision, and the thermally averaged cross sections for s-wave production with Hulthen potential and bound state formation (BSF) for the states 1s1, 1s32 1s5, 2s1, 2s3, 2s5, 2p1, 2p3, 2p5

Basic_definitions.py - contains constant definitions and import Bessel functions

Quintuplet_definitions.py - contains the Quintuplet_DM class definition

BoundStates_definitions.py - contains the Quintuplet_BS class definition

NDE_solver.py - contains the ODE solver function, and does the mass scan for different cases (IMPLEMENTED: Hulthen_free, Hulthen_1s_effective) 
