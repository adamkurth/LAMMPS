#!/bin/sh
#SBATCH --job-name=sim_slurm_tag
#SBATCH --output=./sim_slurm_tag/%x.out
#SBATCH --error=./sim_slurm_tag/%x.err
#SBATCH --time=4
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=
#SBATCH --partition=htc
#SBATCH --qos=normal
#SBATCH --chdir=/home/amkurth/Development/lammps_project/LAMMPS/tutorials/graphene_sheet
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

lmp -in input_cnt.lammps > ./sim_slurm_tag/sim_slurm_tag.out 2> ./sim_slurm_tag/sim_slurm_tag.err
