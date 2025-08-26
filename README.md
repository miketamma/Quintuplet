# Computation of the Quintuplet relic thermal mass

Precise calculation of the relic mass of the Dark Matter, assumed as being part of the 5 representation of the Standard Model SU(2) gauge group.

If you use this code or its results in your work please cite: ...

Useful references: "*Cosmological Implications of Dark Matter Bound States*", 1702.01141

Packages needed:

- SciPy
- NumPy
- mpmath
- time
- matplotlib (plotting not implemented yet, so not needed as of now)


## TL;DR

The final result can be obtained by simply running the 'Boltzmann_Solver.py' script; it will ask as input the number of mass points to scan in the M<sub>DM</sub> ∈ [3,15] TeV interval, and which case to run (more details below). It outputs the list [M<sub>DM</sub>, Ω<sub>DM</sub> h<sup>2</sup>] as a .txt file and optionally it prints it on terminal.

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
- **Boltzmann_Solver.py** (script)
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


### Basic_definitions.py 

Contains constant definitions, imports Bessel functions and d.o.f. tables and defines the basic function to compute Ω<sub>DM</sub> h<sup>2</sup>. 

It also defines the **NDE_solver** function, which gets as inputs the Boltzmann equations and the number of bound states included, and integrates the equations in the interval z ∈ [3, 10<sup>7</sup>], using the SciPy function 'solve_ivp' with the 'Radau' method. Initial values of the DM and BS yields are set to Y = 10<sup>-30</sup>.

### Quintuplet_definitions.py

Defines gauge boson mass with thermal corrections, and imports the thermally averaged DM cross sections

Mainly, it defines the **Quintuplet_DM** class, which takes as input parameter only the DM mass, and initializes a series of functions relevant for the DM.


### BoundStates_definitions.py 

Similarly, it imports and defines bound states functions.

The **Quintuplet_BS** inherits the **Quintuplet_DM** class, and defines additional functions for the BS properties. It takes the following list of inputs:

[ M<sub>DM</sub>, g<sub>I</sub>, n, ℓ, I, BSF function, Group factor constant ]

### Free_Boltz_eqs.py

*Add description*

### Hulthen_eff_Boltz_eqs.py

*Add description*

### Hulthen_network_Boltz_eqs.py

*Add description*

### Boltzmann_Solver.py

This scripts takes as input all the previous ones and runs the required computations.

First, it asks as input the number of mass to scan. The mass range is fixed to be [3, 15] TeV (*might change in the future, not a priority now*).

Second, it asks for the case to run. The implemented cases are:

    1. Tree level (No Sommerfeld enhancement) (~0.6 s/mass)

    2. Sommerfeld enhancement with Hulthen potential (~1.0 s/mass)

    3. Effective Boltzmann equation with 1s bound states (~1.8 s/mass)

    4. Effective Boltzmann equation with 1s, 2s, 2p bound states (~3.5 s/mass)

    5. Network of Boltzmann equations with 1s1 bound states (~3.8 s/mass)

    6. Network of Boltzmann equations with 1s bound states (~8.6 s/mass)

    7. Network of Boltzmann equations with 1s, 2s, 2p bound states (~30 s/mass)

    0. Solve all of them (could take some time)

Finally, it asks if to write the results on terminal or not.

### Results

Output files are saved in this folder (currently only .txt files with the mass scans). 

## ToDo List

- Complete the current README file

- Add possibility to decide the range of masses to scan

- Add possibility to extract directly the thermal mass from the mass scan

- Change output format file? Currently using .txt format to have a simple list to import in Mathematica and make plots

- Implement the possibility to plot results directly from the python script
