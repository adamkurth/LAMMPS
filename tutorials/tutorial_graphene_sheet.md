
### Graphene Sheet and CNT
 
*Follow along* [here](https://lammpstutorials.github.io/lammpstutorials-version1.0/tutorials/04-Graphene.html)

**Goal** to use molecular dynamics and simulate carbon-based structures, and to study the mechanical properties of a graphene sheet. There are two parts to this tutorial.

1. *Graphene Sheet*: 
    - The first part is to generate a graphene sheet using VMD-topotool and then deformed using an applied displacement (method called *out-of-equilibrium  molecular dynamics*)
2. *CNT*: 
   - The second part is to create a single CNT. In case, a reactive force field will be used, and the breaking of bonds under extreme deformation will be simulated.

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

9. Save the `carbon.data` file in the same folder as your future LAMMPS input script. We ar eedont with system generation, and can start on the molecular dynamics simulation.

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

