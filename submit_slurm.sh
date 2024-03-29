#!/bin/bash                                                                     

# Check if the correct number of arguments is provided
if [ "$#" -ne 7 ]; then                                                         
    echo "Usage: $0 <simulation_name> <input_script> <tasks> <partition> <qos> <time> <tag>"             
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

echo "Loading LAMMPS module..."
module load lammps/29Sep2021                                                    
                                                                                
echo "Changing to LAMMPS project directory..."
cd /home/amkurth/Development/lammps_project/LAMMPS/tutorials/graphene_sheet    

# Define simulation parameters                                                  
RUN="sim_${NAME}${TAG}"  # Append TAG to RUN if provided                        
OUTPUT_DIR="./${RUN}"                                                           
echo "Simulation Name: $NAME"
echo "Input Script: $INPUT_SCRIPT"
echo "Tasks: $TASKS"
echo "Partition: $PARTITION"
echo "QOS: $QOS"
echo "Time: $TIME"
echo "Tag: $TAG"
echo "Output Directory: $OUTPUT_DIR"

# Create output directory                                                       
mkdir -p "$OUTPUT_DIR"                                                          

# Create the SLURM job script                                                   
SLURMFILE="${RUN}.sh"                                                           
echo "#!/bin/sh" > $SLURMFILE                                                   
echo "#SBATCH --job-name=$RUN" >> $SLURMFILE                                    
echo "#SBATCH --output=$OUTPUT_DIR/%x.out" >> $SLURMFILE                        
echo "#SBATCH --error=$OUTPUT_DIR/%x.err" >> $SLURMFILE                         
echo "#SBATCH --time=$TIME" >> $SLURMFILE                                       
echo "#SBATCH --ntasks=$TASKS" >> $SLURMFILE                                    
echo "#SBATCH --cpus-per-task=$SLURM_CPUS_PER_TASK" >> $SLURMFILE  # Assumes you want to set CPUs per task to match TASKS
echo "#SBATCH --partition=$PARTITION" >> $SLURMFILE                             
echo "#SBATCH --qos=$QOS" >> $SLURMFILE                                         
echo "#SBATCH --chdir=$PWD" >> $SLURMFILE                                       
echo "export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK" >> $SLURMFILE  # Set OMP_NUM_THREADS to match SLURM_CPUS_PER_TASK
echo "" >> $SLURMFILE                                                           

# Ensure the LAMMPS executable path is correct for your environment             
LAMMPS_EXEC="lmp"  # Use the correct executable name as per your module load

# Add the command to run LAMMPS                                                 
echo "$LAMMPS_EXEC -in $INPUT_SCRIPT > $OUTPUT_DIR/$RUN.out 2> $OUTPUT_DIR/$RUN.err" >> $SLURMFILE

# Submit the job                                                                
echo "Submitting SLURM job: $RUN with $TASKS tasks..."
sbatch $SLURMFILE                                                               
                                                                                
echo "Submitted SLURM job: $RUN"                                                
