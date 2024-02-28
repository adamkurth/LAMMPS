o "changing to LAMMPS directory ------------------------------"
cd /home/amkurth/Development/lammps_project

if [ "$#" -ne 5 ]; then
    echo "Usage: $0 NAME TASKS PARTITION QOS TIME"
    exit 1
fi

# Extract arguments
NAME="$1"
TASKS="$2"
PARTITION="$3"
QOS="$4"
TIME="$5"

# Define simulation parameters
RUN="${NAME}"
SLURMFILE="${RUN}/${RUN}.sh"
OUTPUT_DIR="./${RUN}"
INPUT_SCRIPT="in.graphene"  # Directly specifying the input script name

# Create output directory
if [ -d "$OUTPUT_DIR" ]; then
    echo "Output directory $OUTPUT_DIR already exists. Exiting."
    exit 1
else
    mkdir -p "$OUTPUT_DIR"
fi

# Create the SLURM job script with job parameters
echo "#!/bin/sh" > $SLURMFILE
echo "#SBATCH --job-name=$RUN" >> $SLURMFILE
echo "#SBATCH --output=$OUTPUT_DIR/%x.out" >> $SLURMFILE
echo "#SBATCH --error=$OUTPUT_DIR/%x.err" >> $SLURMFILE
echo "#SBATCH --time=$TIME" >> $SLURMFILE
echo "#SBATCH --ntasks=$TASKS" >> $SLURMFILE
echo "#SBATCH --partition=$PARTITION" >> $SLURMFILE
echo "#SBATCH --qos=$QOS" >> $SLURMFILE
echo "#SBATCH --chdir=$OUTPUT_DIR" >> $SLURMFILE
echo "" >> $SLURMFILE

# Ensure the LAMMPS executable path is correct for your environment
LAMMPS_EXEC="lmp"  # Adjust this to match your actual LAMMPS executable

# Add the command to run LAMMPS, ensuring output and log files are managed within the specified directory
echo "$LAMMPS_EXEC -in $INPUT_SCRIPT -log log.lammps > $RUN.out 2> $RUN.err" >> $SLURMFILE

# Submit the job
sbatch $SLURMFILE

echo "Submitted SLURM job: $RUN ------------------------------"

