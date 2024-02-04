#!/bin/bash

# For running a LAMMPS simulation of a graphene sheet and carbon nanotube

# Load LAMMPS module or set environment variables if needed
# module load lammps  # Uncomment and adjust if you use a module system
# or export PATH=/path/to/lammps:$PATH  # Adjust to your LAMMPS installation

# Assume LAMMPS and all necessary files are in the current directory

# For submitting a LAMMPS job to a queueing system:
# module load lammps  # Uncomment and adjust if you use a module system

# Set the LAMMPS executable name according to your installation
LAMMPS_EXEC=lmp_serial  # Use lmp_mpi for MPI execution if applicable

# Set the filename of your LAMMPS input script
INPUT_SCRIPT=input_cnt.lammps

# Define the log file name
OUTPUT_LOG=cnt_simulation.log

echo "Starting LAMMPS simulation for Carbon Nanotube..."

# Run the LAMMPS simulation
$LAMMPS_EXEC -in $INPUT_SCRIPT -log $OUTPUT_LOG

# Check if the simulation completed successfully
if [ $? -eq 0 ]; then
    echo "LAMMPS simulation completed successfully."
    # Optional: Add any post-simulation commands here, e.g., data analysis
else
    echo "Error: LAMMPS simulation did not complete successfully. Check $OUTPUT_LOG for details."
fi

echo "Script execution finished."

# Optional: Post-processing or analysis commands
# e.g., python analysis_script.py  # Adjust to your post-processing script

echo "Script finished."
