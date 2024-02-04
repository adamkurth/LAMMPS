#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <simulation_name> <input_script>"
    exit 1
fi

# Extract arguments
SIM_NAME="$1"
INPUT_SCRIPT="$2"

# Define simulation parameters
RUN="sim_${SIM_NAME}"
OUTPUT_DIR="./${RUN}"
LAMMPS_EXEC="lmp_serial"  # Adjust if your LAMMPS executable name is different

echo "Simulation Name: $SIM_NAME"
echo "Input Script: $INPUT_SCRIPT"
echo "Output Directory: $OUTPUT_DIR"

# Create a new output directory for this run
mkdir -p "$OUTPUT_DIR"

# Define the log file name within the output directory
OUTPUT_LOG="$OUTPUT_DIR/${SIM_NAME}_simulation.log"
echo "Log file will be written to: $OUTPUT_LOG"

# Construct the LAMMPS command
command="$LAMMPS_EXEC -in $INPUT_SCRIPT -log $OUTPUT_LOG"

# Execution
echo "Starting LAMMPS simulation for $SIM_NAME..."
eval $command
# $LAMMPS_EXEC -in $INPUT_SCRIPT -log $OUTPUT_LOG

# Check if the simulation completed successfully
if [ $? -eq 0 ]; then
    echo "LAMMPS simulation completed successfully for $SIM_NAME."
else
    echo "Error: LAMMPS simulation did not complete successfully for $SIM_NAME."
    echo "Check the error details in the output log: $OUTPUT_LOG"
fi

# Optional: Here you can add commands to move additional generated files into the output directory
# For example:
# mv ${SIM_NAME}_output.* $OUTPUT_DIR/

echo "LAMMPS simulation script execution finished for $SIM_NAME."
echo "All generated files are located in: $OUTPUT_DIR"

# Optional: Add any post-simulation commands here, e.g., data analysis

echo "All processes for $SIM_NAME completed."
