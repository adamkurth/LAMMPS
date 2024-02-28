#!/bin/bash

# Load LAMMPS module
echo "Loading LAMMPS: lammps/29Sep2021 --------------------------"
module load lammps/unstable-5486896

# Change directory to LAMMPS project
echo "Changing to LAMMPS project directory ------------------------------"
cd /home/amkurth/Development/lammps_project

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 NAME INPUT_SCRIPT"
    exit 1
fi

# Extract arguments
NAME="$1"
INPUT_SCRIPT="$2"

# Define simulation parameters
RUN="${NAME}"
OUTPUT_DIR="./${RUN}"

# Print input arguments
echo "Parameters: NAME=$NAME, INPUT_SCRIPT=$INPUT_SCRIPT"
echo "Running LAMMPS simulation locally ------------------------------"

# Ensure the LAMMPS executable path is correct
LAMMPS_EXEC="lmp"  # Adjust based on your installation

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run LAMMPS simulation
$LAMMPS_EXEC -in $INPUT_SCRIPT > $OUTPUT_DIR/$RUN.out 2> $OUTPUT_DIR/$RUN.err

echo "LAMMPS simulation completed: $RUN ------------------------------"

