### Graphene Sheet and CNT
 
*Follow along* [here](https://lammpstutorials.github.io/lammpstutorials-version1.0/tutorials/04-Graphene.html)

**Goal** to use molecular dynamics and simulate carbon-based structures, and to study the mechanical properties of a graphene sheet. There are two parts to this tutorial.

1. *Graphene Sheet*: 
    - The first part is to generate a graphene sheet using VMD-topotool and then deformed using an applied displacement (method called *out-of-equilibrium  molecular dynamics*)
2. *CNT*: 
   - The second part is to create a single CNT. In case, a reactive force field will be used, and the breaking of bonds under extreme deformation will be simulated.

*If following along, all of the contents here will be in the `LAMMPS/tutorials/graphene_sheet` directory.*

#### Graphene Sheet

To generate the initial atom positions, bonds, and angles, we use `VMD`.

1. Open `VMD` and click `Extensions` -> `Modeling` -> `Nanotube Builder` -> `Graphene Sheet`
2. In the `Graphene Sheet` section, we can change the `Size` in the `x` and `y` direction to follow the tutorial to be 4 and 4 for both x and y. 
3. Click `Generate Sheet(s)`.
4. To change the `Drawing Method` in `Graphics` from `Lines` to `DynamicBonds`. This is optional.
5. At this point this is not a molecular dynamics simulation, but a cloud of dots that look like graphene.
6. To create the bonds and angles, we need to use the `Topotool` plugin in `VMD`, this will estimate which dots must be connected by bonds using a distance criterion.
7. In the `VMD` terminal, set the box dimensions bu typing the following:

```bash
molinfo top set a 80  
molinfo top set b 80            
molinfo top set c 80 
```

8. Then generate the LAMMPS file, enter the following:

```bash
topo writelammpsdata carbon.data full
```

- Note: the keyword here `full` corresponds to the LAMMPS atom full style (other possibilites are `atomic`, `bond`, and `charge`). The parameters of the contstraints (bond, length, dihedral, etc.) are also written in the file.
- I moved this to `LAMMPS/tutorials/graphene_sheet/` under the name `carbon.data`.
  
- This file contains information about the positions of the carbon atoms, and the identity of the atoms linked by bonds, angles, dihedrals, and impropers contraints. 

9. Save the `carbon.data` file in the same folder as your future LAMMPS input script. We are done with system generation, and can start on the molecular dynamics simulation.

- Note: I will copy this file and put it in the `LAMMPS/tutorials/` directory as well.

#### LAMMPS Input Script

1. Create a new `.txt` file and name this `input.lammps`. Copy the following lines into the file:

```bash
# Initialisation

variable T equal 300

units real
atom_style full
boundary p p p
pair_style lj/cut 14

bond_style harmonic
angle_style harmonic
dihedral_style opls
improper_style harmonic

special_bonds lj 0.0 0.0 0.5

# System definition
read_data carbon.data
```

Notes: 
- Most of these command lines have been seen in previous tutorials, but there are a few differences. 
- First, the pair style here is `lj/cut` with parameter 14, which is a Lennard-Jones potential with a cutoff. This means that the atoms closer than 14 Angstroms from eachother interact through a Lennard-Jones potential. Notice there are no Columbic interactions, because all the atoms in the pure graphene have charge of 0. The bond, angle, dihedral, and improper styles specify the different potentials used to restrain the positions of the atoms.

- *Interaction between neighbors atoms*: Atoms connected by bond do not typically interact through Lennard-Jones interaction. This ensures here by the `special_bonds` comannd. The three neighbors of the `special_bonds` command are weight factors for the Lennard-Jones interaction between atoms connected by bond (respectively directly bounded, seperated by two bonds, etc.). Example, the first weight factor, with a value of 0, imposes that the two atoms connected by a bond do not interact through Lennard-Jones potential (thus, only interact through bonded potentials).

- the last command, `read_data` imports the `carbon.data` file previously generated with `VMD` which contains the information about the box size, atom positions, etc.

##### Parameters for graphene

Next, the parameters  of both the bonded and non-bonded interactions must be specified. Create a new text file called `PARM.lammps` and copy the following lines into the file:

```bash
pair_coeff 1 1 0.066047 3.4
bond_coeff 1 469 1.4
angle_coeff 1 63 120
dihedral_coeff 1 0 7.25 0 0
improper_coeff 1 5 180
```

- The `pair_coeff` command sets the Lennard-Jones parameters ($ɛ, σ$) for the only type of atoms of the simulation: carbon atom type 1.
- The `bond_coeff` provides the equilibrium distance, $r_0$ and the spring constant $K$ for the harmonic potential imposed between two neighboring carbon atoms, where the potential is :

    $$ E = K(r-r_0)^2 $$

- The `angle_coeff` command gives the equilibrium angle $\theta_0$ and constant for the harmonic potential between three neighboring carbon atoms. The potential is:

    $$ E = K(\theta-\theta_0)^2 $$

- The `dihedral_coeff` and `improper_coeff` give the potential for the constraints between 4 atoms. The `PARAM.file` can be included in the simulation by adding the following line to the `input.lammps` file:

```bash
include	PARM.lammps
```

#### Seperate the atoms into 3 groups

<!-- 480 342 1735 -->
FAX: 480 301 6780 

The goal of the present simulation is to impose a deformation to the sheet. To do this, we isolate the atoms of the two edges of the graphene sheets into two groups, and the displacement will be applied to the atoms of the edge. Add the following lines to the `input.lammps` file:

```bash
# Simulation settings

group gcar type 1
variable xmax equal bound(gcar,xmax)-0.5
variable xmin equal bound(gcar,xmin)+0.5
region rtop block ${xmax} INF INF INF INF INF
region rbot block INF ${xmin} INF INF INF INF
region rmid block ${xmin} ${xmax} INF INF INF INF
```

The first command includes all of the atoms of type 1 (i.e. all that atoms here) in a group named `gcar`. Then, two variables are defined:  

$x_{max}$ corresponds to the coordinate of the last atoms along the $x$ minus 0.5 Angstroms, and $x_{min}$ to the coordinate of the firsst atom along the $x$ axis plus 0.5 Angstroms. Then 3 regions are defined, and correspond respectively to: 
$$
\begin{align*}
    x &< x_{min} \\
    x_{min} > x &> x_{max} \\
    x &> x_{max}
\end{align*}
$$


Then let us define 3 groups of atoms corresponding to the atoms located in each of the 3 regions respectively. Add the following lines to the `input.lammps` file:

```bash
group gtop region rtop
group gbot region rbot
group gmid region rmid
```

Note: when running the simulation, the number of atoms in each group is printed in the terminal, (and log.lammps file), it's a way of controlling that no mistake was made in the definition of the groups. e.g. `20 atoms in group gtop`.

#### Thermalisation and dynamics

Add the following lines to the `input.lammps` file:

```bash
velocity gmid create ${T} 48455 mom yes rot yes
fix mynve all nve
compute Tmid gmid temp
fix myber gmid temp/berendsen ${T} ${T} 100
fix_modify myber temp Tmid
```

- The `velocity create` command gives the initial velocities of the atoms of the group gmid, assuring an initial temperature of 300 K. The `mom` and `rot` keywords are used to specify that the linear and angular momentum of the group gmid are conserved. 
- NVE fix is appied to all atoms, thus ensuring that atoms position are recalculated in time, and Brentsen thermostat is applied to the group gmid, to ensure that the temperature of the group gmid is maintained at 300 K.
- The `fix_modify` ensures that the fix Berendsen uses the temperature of the group gmid as an input, instead of the temperatrue of the whole system (global temp). The atoms of the edges are not thermalised because their motion will be restrained in the next part of the input script.

#### Restrain the motion of the edges

To restrain the motion of the atoms at the edges, add the following lines to the `input.lammps` file: 

```bash
fix mysf1 gtop setforce 0 NULL 0
fix mysf2 gbot setforce 0 NULL 0
velocity gtop set 0 NULL 0
velocity gbot set 0 NULL 0
```

- the two `setforce` commands cancel the forces applied on the atoms of the two edges, respectively, during the whole simulation along the \(x\) and \(z\), and the velocity commands set the initial velocity along the \(x\) and \(z\) to 0 for the atoms of the edges. Thus, the atoms of the edges will remain immobile during the simulation (or they would if no other command was added).


#### Data Extraction

In order the measure the straing and stress of the graphene sheet, lets extract the distance \(L\) between the two edges as well as the force applied on the edges. Let's also add a command to print the atoms coordinates in the `lammpstrj` file, for every 1000 timesteps. Add the following lines to the `input.lammps` file:

```bash
variable	L equal xcm(gtop,x)-xcm(gbot,x)
fix 		at1 all ave/time 10 100 1000 v_L file length.dat
fix 		at2 all ave/time 10 100 1000 f_mysf1[1] f_mysf2[1] file force.dat
dump 		mydmp all atom 1000 dump.lammpstrj
```

- Note that the values of the force on each step are extracted from the fixes `setforce` `mysf1` and `mysf2` and are stored in the file `force.dat`. By calling them `f_` the same way variables are called using `v_` and computes are called using `c_`. 
- A fix `setforce` cancels all the forces on a group of atoms every timestep, but allows one to extract the values of the force before cancelation. 
- The `ave/time` command computes the average of the force on the atoms of the edges every 1000 timesteps, and stores the result in the file `force.dat`
- The `xcm` command computes the center of mass of the group `gtop` and `gbot` along the x axis, and the difference between the two is stored in the file `length.dat`. The `dump` command writes the coordinates of the atoms in the file `dump.lammpstrj` every 1000 timesteps.
- The `dump` command writes the coordinates of the atoms in the file `dump.lammpstrj` every 1000 timesteps.


#### Run

With one small line added to add an equilibrium step, the input script is now complete. Add the following lines to the `input.lammps` file:

```bash
# eqilibrium
thermo 100
thermo_modify temp Tmid

# Run

timestep 1.0
run 5000
```

with the `thermo_modify` command, we specify to LAMMPS that we want the entire temperature \( T_{mid}\) to be printed in the terminal, not the temperature of the entire system (because of frozen edges, the temp of the entire system is not relavant). Then let us perform a loop. At each step of the loop, the edges are slightly displaced, and the simulation runs for a short time. 

```bash
variable var loop 10
label loop
displace_atoms gtop move 0.1 0 0
displace_atoms gbot move -0.1 0 0
run 1000
next var
jump input.lammps loop
```

- What you observe should resemble the following: [Click here](https://www.youtube.com/embed/o5IoCVWpPKg)
- The sheet is progressively elongated, and the carbon honeycombs are being deformed. You can increase the number of iterations of the loop (variable `var`) to force a longer elongation.

#### *Notes*

- Always remember that what you measure and observe is only as good as your force field. 
- With the present force field, no matter how large is the imposed deformation, the bonds will never break. To study such bond breaking, one has to use a reactive force field, sucha as in the following step. 

___


#### Carbon Nanotube

Using the same protocol as for the graphene sheet, generate a carbon nanotube datafile, call this `cnt.data`. Before generating CNT, untick the `bonds`.

1. To do this hit `Extensions` -> `Modeling` -> `Nanotube Builder`. Then, in the `Nanotube Builder` section, untick `bonds` and enter the desired parameters for the CNT.

2. Then create a LAMMPS input file, and type the following:
``` bash
touch cnt.data
```

- Inside of `LAMMPS/tutorials/graphene_sheet` directory. 

3. Then, copy the following lines into the new input: 

```bash
# Initialisation

variable T equal 300

units metal
atom_style full
boundary p p p
pair_style airebo 2.5 1 1

# System definition
read_data cnt.data
```

- Can rename the name in the `VMD` window to `CNT`.

Differences: 

- the first difference with the previous case (graphene) is the units `metal`, instead of `real`. This is a choice that is imposed by the force field we are going to use (careful, the time is in pico seconds \( 1e -12 s\)) with `metal` instead of femto seconds \(1e-15s\) with `real`.

- The second difference is the `pair_style` command, which is `airebo` is a reactive force field. With reactive force field, bonds between atoms are deduced in real time according to the distance between atoms. 

4. Import the LAMMPS data file and set the `pair_coeff` command.
    
```bash
read_data carbon.data
pair_coeff * * CH.airebo C
```

- The `CH.airebo` file contains the parameters of the reactive force field. And can be downloaded [here](https://github.com/lammps/lammps/blob/develop/potentials/CH.airebo). 
- All the rest of the script is very similar to the previous one, except that the `fix` commands are not used to restrain the motion of the edges, because the reactive force field will take care of the bond breaking.

5. Enter the following lines after downloading the `CH.airebo` file:

```bash
# Simulation settings

group gcar type 1
variable zmax equal bound(gcar,zmax)-0.5
variable zmin equal bound(gcar,zmin)+0.5
region rtop block INF INF INF INF ${zmax} INF
region rbot block INF INF INF INF INF ${zmin}
region rmid block INF INF INF INF ${zmin} ${zmax}

group gtop region rtop
group gbot region rbot
group gmid region rmid

velocity gmid create ${T} 48455 mom yes rot yes
fix mynve all nve
compute Tmid gmid temp
fix myber gmid temp/berendsen ${T} ${T} 0.1
fix_modify myber temp Tmid
```

6. For a change, lets impose a constant velocity to the atoms of one edge, while maintaining the other edge fix. To do so, one needs to cancel the forces (thus the acceleration) on the atoms of the edges using the `setforce` command, and set the value of the velocity along the \(z\) direction.

First, as an equlibrium step, lets set the velocity to 0. Add the following lines to the `input_cnt.lammps` file:

```bash
fix mysf1 gbot setforce NULL NULL 0
fix mysf2 gtop setforce NULL NULL 0
velocity gbot set NULL NULL 0
velocity gtop set NULL NULL 0

variable pos equal xcm(gtop,z)
fix at1 all ave/time 10 100 1000 v_pos file cnt_deflection.dat
fix at2 all ave/time 10 100 1000 f_mysf1[1] f_mysf2[1] file force.dat
dump mydmp all atom 1000 dump.lammpstrj

thermo 100
thermo_modify temp Tmid

# Run
timestep 0.0005
run 5000
```

- At the start of the equilibrium step, you can see that the temperature deviates from the target temperature of 300 K (after a few picoseconds the temperature reaches the target value): 

- Output:

```bash
Step          Temp          E_pair         E_mol          TotEng         Press     
0   300           -5084.7276      0             -5058.3973     -1515.7017    
100   237.49462     -5075.4114      0             -5054.5671     -155.05545    
200   238.86589     -5071.9168      0             -5050.9521     -498.15029    
300   220.04074     -5067.1113      0             -5047.7989     -1514.8516    
400   269.23434     -5069.6565      0             -5046.0264     -174.31158    
500   274.92241     -5068.5989      0             -5044.4696     -381.28758    
600   261.91841     -5065.985       0             -5042.9971     -1507.5577    
700   288.47709     -5067.7301      0             -5042.4111     -312.16669    
800   289.85177     -5066.5482      0             -5041.1086     -259.84893    
900   279.34891     -5065.0216      0             -5040.5038     -1390.8508    
1000   312.27343     -5067.6245      0             -5040.217      -465.74352
```

7. Change the velocity to 30 m/s and run for longer time by add the following lines to the `input_cnt.lammps` file: 

- Add this to the very bottom of the input script.

```bash
# 0.15 A/ps = 30 m/s
velocity gtop set NULL NULL 0.15
run 250000
```

