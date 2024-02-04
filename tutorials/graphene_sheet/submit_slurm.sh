#!/bin/bash                                                                     

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <simulation_name> <input_script>"
    exit 1
fi

# Extract arguments
SIM_NAME="$1"
INPUT_SCRIPT="$2"

# Load LAMMPS module
echo "Loading LAMMPS module..."
module load lammps/29Sep2021

# Change directory
echo "Changing to the LAMMPS project directory..."
cd /home/amkurth/Development/lammps_project/LAMMPS/tutorials/graphene_sheet

# Define simulation parameters
RUN="sim_${SIM_NAME}"
OUTPUT_DIR="./${RUN}"
echo "Simulation Name: $SIM_NAME"
echo "Input Script: $INPUT_SCRIPT"
echo "Output Directory: $OUTPUT_DIR"

# Create a new output directory for this run
mkdir -p "$OUTPUT_DIR"

# Define SLURM job script name
SLURMFILE="${RUN}.sh"

# Create the SLURM job script
echo "#!/bin/sh" > $SLURMFILE
echo "#SBATCH --job-name=$RUN" >> $SLURMFILE
echo "#SBATCH --output=$OUTPUT_DIR/%x.out" >> $SLURMFILE
echo "#SBATCH --error=$OUTPUT_DIR/%x.err" >> $SLURMFILE
# Adjust the following SLURM parameters as necessary for your system and requirements
echo "#SBATCH --time=01:00:00" >> $SLURMFILE
echo "#SBATCH --ntasks=1" >> $SLURMFILE
echo "#SBATCH --partition=standard" >> $SLURMFILE
echo "#SBATCH --qos=standard" >> $SLURMFILE
echo "#SBATCH --chdir=$PWD" >> $SLURMFILE
echo "" >> $SLURMFILE

# Add the command to run LAMMPS
echo "lmp_serial -in $INPUT_SCRIPT > $OUTPUT_DIR/$RUN.out 2> $OUTPUT_DIR/$RUN.err" >> $SLURMFILE

# Submit the job
echo "Submitting SLURM job..."
sbatch $SLURMFILE
echo "Submitted SLURM job: $RUN"
