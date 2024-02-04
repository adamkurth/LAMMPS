#!/bin/bash

# Correct the number of expected arguments in the initial check
if [ "$#" -ne 7 ]; then
    echo "Usage: $0 NAME INPUT_SCRIPT TASKS PARTITION QOS TIME TAG"
    exit 1
fi

# Extract arguments
NAME="$1"
INPUT_SCRIPT="$2"
TASKS="$3"
PARTITION="$4"
QOS="$5"
TIME="$6"
TAG="$7"

# Define simulation parameters
RUN="sim_${NAME}${TAG}"  # Append TAG to RUN if provided
SLURMFILE="${RUN}.sh"
OUTPUT_DIR="./${RUN}"

# Create the SLURM job script
echo "#!/bin/sh" > $SLURMFILE
echo "#SBATCH --job-name=$RUN" >> $SLURMFILE
echo "#SBATCH --output=$OUTPUT_DIR/%x.out" >> $SLURMFILE
echo "#SBATCH --error=$OUTPUT_DIR/%x.err" >> $SLURMFILE
echo "#SBATCH --time=$TIME" >> $SLURMFILE
echo "#SBATCH --ntasks=$TASKS" >> $SLURMFILE
echo "#SBATCH --partition=$PARTITION" >> $SLURMFILE
echo "#SBATCH --qos=$QOS" >> $SLURMFILE  # Add QOS to SLURM script
echo "#SBATCH --chdir=$PWD" >> $SLURMFILE
echo "" >> $SLURMFILE

# Ensure the LAMMPS executable path is correct for your environment
LAMMPS_EXEC="lmp_serial"  # Adjust for your LAMMPS installation, e.g., lmp_mpi for MPI

# Add the command to run LAMMPS
echo "$LAMMPS_EXEC -in $INPUT_SCRIPT > $OUTPUT_DIR/$RUN.out 2> $OUTPUT_DIR/$RUN.err" >> $SLURMFILE

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Submit the job
sbatch $SLURMFILE

echo "Submitted SLURM job: $RUN"
