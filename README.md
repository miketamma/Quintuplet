# Computation of the Quintuplet relic thermal mass

Precise calculation of the relic mass of the Dark Matter, assumed as being part of the 5 representation of the Standard Model SU(2) gauge group.


## TL;DR

The final result can be obtained by simply running the 'Boltzmann_solver.py' script; it will ask as input the number of mass points to scan in the M<sub>DM</sub> ∈ [3,15] TeV interval, and which case to run (more details below). It outputs the list [M<sub>DM</sub>, Ω<sub>DM</sub> h<sup>2</sup>] as a .txt file and optionally it prints it on terminal.

## Structure of the code

The code is structured in a subset of files and scripts that contain the different definitions necessary to compute the Boltzmann equations.

These are described below:

- **Num_lists** (folder)
- **Basic_definitions.py** (script)
- **Quintuplet_definitions.py** (script)
- **BoundStates_definitions.py** (script)
- **Free_Boltz_eqs.py** (script)
- **Hulthen_eff_Boltz_eqs.py** (script)
- **Hulthen_network_Boltz_eqs.py** (script)
- **Results** (folder)


### Num_lists

This folder contains all the necessary pretabulated quantities as .txt files. Scripts to (re)compute these tables are also provided.

- 'K1_hiprec.txt' and 'K2_hiprec.txt' contain the tabulated Log10 of the modified Bessel functions of first and second kind, K<sub>1</sub>(z) and K<sub>2</sub>(z), as function of the variable z ∈ [1, 10<sup>9</sup>].  
- 'nDMeq_z.txt' and 'YDMeq_z.txt' contain the tabulated Log10 of the functions K<sub>2</sub>(z)/z and K<sub>2</sub>(z) * z<sup>2</sup>, as function of the variable z ∈ [1, 10<sup>9</sup>].

The above tables are computed by the **Hiprec_Bessels.py** script

- 'gstarsT.csv' and 'sqrtgsT.csv' contain the number of degrees of freedom and the square root of the effective number of d.o.f., as function of the temperature T ∈ [10<sup>-5</sup>, 10<sup>5</sup>] GeV.

- 'NoSommTA.txt' contains the tabulated Log10 of the thermal average of the Dark Matter annihilation cross section for isospin states I = 1, 3, 5, without any Sommerfeld enhancement, as function of z ∈ [10<sup>0</sup>, 10<sup>6</sup>] and epsilon T ∈ [10<sup>-5</sup>, 10<sup>2</sup>]. 

- 'SommHueltenTA.txt' contains the tabulated Log10 of the thermal average of the Dark Matter annihilation cross section for isospin states I = 1, 3, 5, including the Sommerfeld enhancement from a Hulthen potential, as function of z ∈ [10<sup>0</sup>, 10<sup>6</sup>] and epsilon T ∈ [10<sup>-5</sup>, 10<sup>2</sup>].

The thermal averages above are computed by the **Thermal_Averages.py** script. Note that the final result is normalized by the 'reference cross section' σ<sub>0</sub>.

- 'BSF_xxx_TA.py' contains the tabulated Log10 of the thermal average of the bound state formation cross sections, as function of z ∈ [10<sup>0</sup>, 10<sup>6</sup>] and epsilon T ∈ [10<sup>-5</sup>, 10<sup>2</sup>]. The bound states are xxx = 1s1, 1s3, 1s5, 2s1, 2s3, 2s5, 2p1, 2p3, 2p5.

The thermal averages above are computed by the **Thermal_Averages_bsf.py** script. Note that the final result is normalized by the 'reference cross section' σ<sub>0</sub><sup>'</sup>.

In Num_lists folder there are the txt files for K1 and K2 (Bessel functions), computed with high precision, and the thermally averaged cross sections for s-wave production with Hulthen potential and bound state formation (BSF) for the states 1s1, 1s32 1s5, 2s1, 2s3, 2s5, 2p1, 2p3, 2p5

Basic_definitions.py - contains constant definitions and import Bessel functions

Quintuplet_definitions.py - contains the Quintuplet_DM class definition

BoundStates_definitions.py - contains the Quintuplet_BS class definition

NDE_solver.py - contains the ODE solver function, and does the mass scan for different cases (IMPLEMENTED: Hulthen_free, Hulthen_1s_effective) 
