#!/bin/bash

# Check if sufficient arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 NAME INPUT_SCRIPT"
    exit 1
fi

# Extract arguments
NAME="$1"
INPUT_SCRIPT="$2"

# Define simulation parameters
RUN="sim_${NAME}"
OUTPUT_DIR="./${RUN}"
LAMMPS_EXEC="lmp_serial"  # Adjust if your LAMMPS executable name is different

# Create a new output directory for this run
mkdir -p "$OUTPUT_DIR"

# Construct the LAMMPS command
command="$LAMMPS_EXEC -in $INPUT_SCRIPT"
command="$command > $OUTPUT_DIR/$RUN.out 2> $OUTPUT_DIR/$RUN.err"

# Execution
echo "Running LAMMPS simulation: $NAME..."
eval $command

echo "LAMMPS simulation completed for $NAME. Output and Error logs are in $OUTPUT_DIR."
