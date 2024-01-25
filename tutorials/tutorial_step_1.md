## LAMMPS Tutorial Step: Input Script

## Part A: Simple Molecular Dynamics

*Followed along here using* [LAMMPS Tutorials](https://lammpstutorials.github.io/lammpstutorials-version1.0/tutorials/01-SimpleMolecularSimulation.html#input)

1. We make a new file called `input_01.lammps` by running the following command:

```bash
touch input_01.lammps
```
2. We open the file in a text editor, and paste the following:

```bash
# PART A - ENERGY MINIMIZATION
# 1) Initialization
# 2) System definition
# 3) Simulation settings
# 4) Visualization
# 5) Run
```

- The `#` symbol is used to comment out lines in the script.

- The script is broken down into 5 parts, and we will go through each part.
  - 1. `Initialization` - This is where we define the units, atom style, and the potential.
  - 2. `System definition` - This is where we define the simulation box, the atoms, and the velocities. 
  - 3. `Simulation settings` - This is where we define the simulation settings, such as the timestep, the number of steps, and the output.
  - 4. `Visualization` - This is where we define the visualization settings, such as the dump file.
  - 5. `Run` - This is where we run the simulation.
  
#### Initialization

By pasting this in the `input_01.lammps` file:

```bash
# 1) Initialization
units lj
dimension 2
atom_style atomic
pair_style lj/cut 2.5
boundary p p p
```

- The first line tells us that we want to use the unit `lj` (Lennard-Jones), for which all quantities are unitless. 
- The second line tells us that the simulation is bidimensional (2D).
- The third line tells us that the atomic style to be used, which means that each atom is represented by a single point in space.
- The fourth line tells that the atoms are going to interact through a Lennard-Jones potential, with a cut-off to 2.5 (unitless). 
- final line tells us that the periodic boundary conditions will be used in all directions. `p p p` means periodic in x, y, and z directions.

Running this command does nothing, but it is a good way to check if the script is correct.

```bash
lmp_serial -in input_01.lammps
```
With an output of:

```bash
LAMMPS (2 Aug 2023)
Total wall time: 0:00:00
```
- A mistake in the script will result in something like this: 

    ```bash
    LAMMPS (20 Nov 2019)
    ERROR: Unknown command: atom_stile	atomic (src/input.cpp:232)
    Last command: atom_stile	atomic
    ```

#### System definition

By pasting this in the `input_01.lammps` file:

```bash
# 2) System definition
region myreg block -30 30 -30 30 -0.5 0.5
create_box 2 myreg
create_atoms 1 random 1500 341341 myreg
create_atoms 2 random 100 127569 myreg
```

- First line creates a region of space, called `myreg` that is a block (rectangular cuboid) of size -30 to 30 along x and y, and -0.5 to 0.5 along z. Note that this is expressed in non-dimensional form due to the `units lj` command. 
- The second line creates a simulation box, called `myreg`, that contains 2 atom types.
- The third line specifies that 1500 atoms of type 1 must be created randomly in the region `myreg`.
- The fourth line specifies that 100 atoms of type 2 must be created randomly in the region `myreg`.
- The integer `341341` is the random seed that can be changed in order to obtain different random distributions of atoms, also for reproducibility purposes.

The output of rerunning this script is:
    
```bash
LAMMPS (20 Nov 2019)
Created orthogonal box = (-30 -30 -0.5) to (30 30 0.5)
1 by 1 by 1 MPI processor grid
Created 1500 atoms
create_atoms CPU = 0.000807692 secs
Created 100 atoms
create_atoms CPU = 4.097e-05 secs
Total wall time: 0:00:00
```

- First, it creates the box, then it creates the atoms. 

#### Simulation settings


By pasting this in the `input_01.lammps` file:

```bash
# 3) Simulation settings
mass 1 1
mass 2 1
pair_coeff 1 1 1.0 1.0
pair_coeff 2 2 0.5 3.0
```
- The two first commands attribute a mass of 1 (unitless) to the two atom types.
- Third line sets the Lennard-Jones parameters coefficients (epsilon and sigma) for the interaction between atoms of type 1.

##### Note: Lennard-Jones parameters

LAMMPS calculates the cross coefficients (of atom types 1,2) using the geometric average (due to lacking rigorous argument): 

- $ ɛ_{ij} = \sqrt{ɛ_{ii}ɛ_{jj}}$, $σ_{ij} = \sqrt{σ_{ii}σ_{jj}}$

- But, $\varepsilon_{ij} = \sqrt{\varepsilon_{ii} \varepsilon_{jj}}$, $\sigma_{ij} = \frac{\sigma_{ii} + \sigma_{jj}}{2}$, are more commonly used. 

- Note: Neither geometric nor the arithmetic rule is base on rigorous argument, so here the geometric will do just fine.

#### Visualization/Run

By pasting this in the `input_01.lammps` file:

```bash
# 4) Visualization
thermo 10

# 5) Run
minimize 1.0e-4 1.0e-6 1000 10000
```

- The `thermo` command tells LAMMPS to print out the thermodynamic properties (temp, energy, pressure, etc.) every 10 timesteps. 
- Second line asks to perform an energy minimization.

##### Note: Energy minimization
Energy minimization procedure consits of adjusting the coordinates of the atoms until one of the stopping criteria is met. Here there are 4 such criteria:
1. Change in energy between two iterations is less than 1.0e-4 (unitless).
2. The maximum force between the atoms in the system are lower than 1.0e-6 (unitless).
3. The maximum number of iterations is 1000.
4. The maximum number of times the force and energy have been evaluated is 10000. 

- The output of the script is:

    ![Alt text](./screenshot_1_24_24.png "Output of Visualization/Run")

- One section of the output looks something like this: 
```bash
Step        Temp    E_pair          E_mol   TotEng      Press
0             0 5.8997404e+14           0 5.8997404e+14 1.5732641e+15
10            0     56275037            0     56275037 1.5007118e+08
20            0    17731.329            0    17731.329    47458.738
30            0    350.68529            0    350.68529    972.14134
40            0    13.745948            0    13.745948    48.748312
50            0    0.5033657            0    0.5033657    8.3398718
60            0   -1.4427524            0   -1.4427524    1.1131474
70            0   -1.7136665            0   -1.7136665 -0.038162464
80            0   -1.7516642            0   -1.7516642  -0.15686171
81            0   -1.7518285            0   -1.7518285  -0.15730928
```
- At step 0, this is the energy of the system, which is very high (5.8997404e+14) (unitless). This is expected due to the random atoms being created at random positions within the simulation box, and some are very close or overlapping with eachother. This result in a large initial energy, repulsive part of the Lennard-Jones interaction potential.
- As the energy minimization proceeds, the energy of the system decreases, and reaches a negiative value, which is more acceptable. This indicates that the atoms have been displaced to more stable and reasonable distances from eachother. 


- 


The timestep is unitless, the temperature is in units of $\varepsilon/k_B$, the energy is in units of $\varepsilon$, and the pressure is in units of $\varepsilon/\sigma^2$.

- The first column is the timestep, the second is the temperature, the third is the potential energy, the fourth is the molecular energy, the fifth is the total energy, and the sixth is the pressure.

This line tells is how many iterations were performed, and the stopping criterion that was met:

```bash
Minimization stats:
  Stopping criterion = energy tolerance
```
This warning is also good to know:

```bash
WARNING: Temperature for thermo pressure is not for group all (src/thermo.cpp:527)
```
Or: 
```bash
WARNING: Using 'neigh_modify every 1 delay 0 check yes' setting during minimization
```

- This is simply an indication that LAMMPS is building the pairwise neighbor lists (that are used by LAMMPS to evaluate the interaction between atoms)

## Part B: Simple Molecular Dynamics

Following here: [Simple molecular dynamics](https://lammpstutorials.github.io/lammpstutorials-version1.0/tutorials/01-SimpleMolecularSimulation.html#input)

Adding to the `input_01.lammps` file:

```bash
# PART B - MOLECULAR DYNAMICS
# 4) Visualization
thermo 1000
variable kinetic_energy equal ke
variable potential_energy equal pe
variable pressure equal press
fix myat1 all ave/time 10 1 10 v_kinetic_energy v_potential_energy v_pressure file energy.dat

# 5) Run
fix mynve all nve
fix mylgv all langevin 1.0 1.0 0.1 1530917
fix myefn all enforce2d
timestep 0.005
run 10000
```
- variables have been defined in order to print out the kinetic and potential energy, and the pressure of the system. This is all outputted in the `energy.dat` file.
- Then the next section, the `fix nve` command is used to update the positions and velocities of the atoms, in the group `all` (*most important command here!*). 
- The second `fix` applies a Langevin thermostat to the atoms in the group of `all`. with the desired temperature of 1. and the damping parameter of 0.1.   
- Number `1530917` is the random seed for the Langevin thermostat.
- Third `fix` command is used to ensure that the atoms remain within the 2D plane.
- Choose the timestep to be 0.005, and ask to run 10000 timesteps. 

Output: 
    
```bash
   Step          Temp          E_pair         E_mol          TotEng         Press     
        81   0             -1.7518285      0             -1.7518285     -0.15730928   
      1000   1.0001038     -1.2883573      0             -0.28887856     1.0156263    
      2000   1.0422003     -1.322024       0             -0.28047505     0.82501505   
      3000   1.0724464     -1.311527       0             -0.23975095     0.79668047   
      4000   1.0102569     -1.3137774      0             -0.30415191     0.76858602   
      5000   1.0200858     -1.3227725      0             -0.30332428     0.73413299   
      6000   0.99159274    -1.2915594      0             -0.30058645     0.90035901   
      7000   0.97611235    -1.3081526      0             -0.3326503      0.78604065   
      8000   1.0088726     -1.3049357      0             -0.29669365     0.76993029   
      9000   0.9959133     -1.3173939      0             -0.32210307     0.76004654   
     10000   1.0152116     -1.3434861      0             -0.32890901     0.62726981   
     10081   1.0203122     -1.3026024      0             -0.28292783     0.75021683   
```
- Second column shoes that the temperature: $T=1$, as requested.
- temperature starts from 0, but reaches the expected value of 1 after a few timesteps.

The output here: 

```bash
Total # of neighbors = 8514
Ave neighs/atom = 5.32125
Neighbor list builds = 1152
Dangerous builds = 0
Total wall time: 0:00:01
```
 
- If `Dangerous builds = 0`, then the simulation is stable. If it is not, then the simulation is unstable, and the timestep should be decreased.

Say this was not the case...

- we should add this line to tell LAMMPS to rebuild the neighbor list every timestep:
  
```bash
# under Simulation settings section...
neigh_modify every 1 delay 5 check yes
```
-  Re-run the simulation, and the output should be:

```bash
Total # of neighbors = 8514
Ave neighs/atom = 5.32125
Neighbor list builds = 1152
Dangerous builds = 0
Total wall time: 0:00:01
```

##### Notes:

- LAMMPS reads top to bottom, so these will be executed after the energy minimization.
- `thermo` command is called a second time, previously with a value of `10`, will be replaced as soon as the new `thermo` command is called.
