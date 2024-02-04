### Submitting Jobs to AGAVE with LAMMPS

Notes for myself to refer to in the future.

1. Login to Agave with ssh and navigate to the directory where the input file is located.

- To find a specific version of lammps, use the following command.
```bash 
module avail # check if lammps is available 
```

- This worked for me.
```bash
module load lammps/29Sep2021 # no errors
cd /home/amkurth/Development/lammps_project
```
- Note that lammps responds to the `lmp` command, instead of `lmp_serial` like on mac.

2. Create a input file for lammps and save it as input.lammps inside of the directory where `run_lammps.sh` is located (should be `lammps_project`). 

```bash
vim input.lammps
```
   - Copy and paste the desired input script into the file, `input.lammps`.

1. Run the following command to submit the job to agave. 

```bash
./run_lammps.sh NAME INPUT_SCRIPT TASKS PARTITION QOS TIME TAG

# Example
./run_lammps.sh lammps_test input.lammps 10 normal htc 4 1
```

