## LAMMPS Tutorial Step 2: Trajectories Visualization

After cloning the github repo here [github repo](https://github.com/lammps/lammps), and following the documentation for Mac OSX from here [docs](https://docs.lammps.org/Install.html). 

## First Simulation: 
### Couette and Poiseuille flow in a 2d channel

*Followed along* [here](https://youtu.be/7RtRerwJqQw?si=_mIPJXG0qJ6AxLLN)

*If following along, all of the contents here will be in the `LAMMPS/tutorials/step_2` directory.*
  
Make sure that the github repo is cloned and ran `cmake` and `make install` in the `build` directory.

1. Go to LAMMPS website and find the movies section. [LAMMPS Movies](https://www.lammps.org/movies.html). 
2. To run `Couette and Poiseuille flow in a 2d channel`, click the link.
3. In the cloned repo, go to `lammps/examples/flow` and find the `in.flow.coutte.lmp` file.
4. Make new directory called `lammps` and `flow` and copy the `in.flow.coutte.lmp` file into it.

- Make sure you're in the `lammps` directory and run the following command:

    ```bash
    # step 4
    cd lammps
    mkdir lammps && cd lammps 
    mkdir flow && cd flow
    cp lammps/examples/flow/in.flow.couette lammps/lammps/flow 
    ```

5. Run this command by the following: 

``` bash
lmp_serial -in in.flow.coutte
```

- Generalizing: You normally run the LAMMPS command in the directory where your input script is located. If you have a script called `in.myinput` in the current directory, you can run it with:

```bash
lmp_serial -in in.myinput
```

7. This will create a `log.lammps` file, where everything from the simulation is stored.

8. To view the simulation, rerun the simulation with the following command:

8. Ensure `vmd` is installed from here, 
    - [VMD](https://www.ks.uiuc.edu/Research/vmd/)
  
9. First, we need to open the `in.flow.coutte` file and change the `dump` command to the following:

    ```bash
    # ... in.flow.coutte
    dump		1 all atom 500 dump.flow #uncomment this line
    ```

10. We remove the commented out line, and if we rerun the simulation, it will create a `dump.flow` file.  

    ``` bash
    lmp_serial -in in.flow.coutte
    ```

- We see that `dump.flow` file has been outputted. 

11. Now open `vmd` and load new molecule, by openoing the `dump.flow` file.

12. Inside `vmd`, click `File` -> `New Molecule` -> find `dump.flow` -> `Open`
13. Now, in the same menu, click the drop down: `Determine file type` -> `LAMMPS Trajectory` -> `Load` 

- Now we can see the simulation.

14. To view the simulation, click `Graphics` -> `Representations` -> `Drawing Method` -> `VDW` -> `Apply` 

- Play around in here, and see what you can do. To view the animation, the play icon is in the bottom right corner.

--- 

## Part 2: Trajectories Visualization

*Followed along here using* [LAMMPS Tutorials](https://lammpstutorials.github.io/lammpstutorials-version1.0/tutorials/01-SimpleMolecularSimulation.html#trajectories)


*For a quick guide using VMD, click* [here](https://lammpstutorials.github.io/lammpstutorials-version1.0/miscellaneous/vmd.html)

To view the simulation, we need to use `vmd` again. And like before, we need to dump the positions of the atoms in a file at a regular interval. Will add from `tutorial_step_1.md`, under the `# visualization` section.

```bash
dump mydmp all atom 1000 dump.lammpstrj
```

Run LAMMPS with the following command:

```bash
# pwd: LAMMPS/tutorials/step_1
lmp_serial -in input_01.lammps
```
- will output a file named `dump.lammpstrj` in the same directory, and can be opened using VMD.

Load this file in VMD, and you will see the simulation.

- To view the simulation, click `Graphics` -> `Representations` -> `Drawing Method` -> `VDW` -> `Apply`. 
