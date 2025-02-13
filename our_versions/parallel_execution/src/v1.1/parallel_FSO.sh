#!/bin/bash

# -------------------------------------------------------------------
# Notes
# Use absolute paths for consistency across nodes
# PBS tells the scheduler how many chunks and cores to allocate for the job
# MPI dictates how many processes to spawn running this code
# PLEASE KEEP IN MIND THAT (PBS_CHUNKS * PBS_CORES_PER_CHUNK >= MPI_RUN_PROCESSES)
# -------------------------------------------------------------------

# max walltime 6h
#PBS -q short_cpuQ

# expected timespan for execution
#PBS -l walltime=06:00:00

# chunks (~nodes) : cores per chunk : shared memory per chunk (?)
#PBS -l select=8:ncpus=5:mem=2gb

readonly C_PROGRAM_PATH=~/hpc_project/our_versions/parallel_execution/src/v1.1/parallel_FSO.c
readonly EXECUTABLE_PATH_AND_NAME=~/hpc_project/our_versions/parallel_execution/src/v1.1/parallel_FSO

# get dependencies
module load mpich-3.2
# build
mpicc $C_PROGRAM_PATH -g -Wall -fopenmp -lm -std=c99 -o $EXECUTABLE_PATH_AND_NAME 
# run

# mpirun.actual -n 5 $EXECUTABLE_PATH_AND_NAME 100 100 100 100
# mpirun.actual -n 1 $EXECUTABLE_PATH_AND_NAME 10 2 100 10


# for i in {1..7}; do
#     for ((i = 1; i <= 20; i=i+1)); do
#         mpirun.actual -n $i $EXECUTABLE_PATH_AND_NAME 100 1000 100 1    
#     done
#         for ((i = 100; i <= 2000; i=i+100)); do
#         mpirun.actual -n 20 $EXECUTABLE_PATH_AND_NAME 20 1000 100 1    
#     done
# done

# # keep the number of schools fixed but increase the number of fishes for each one
# for ((i = 1000; i <= 10000; i=i+1000)); do
#     mpirun.actual -n 20 $EXECUTABLE_PATH_AND_NAME $i 1000 1000 1    
# done

# PLAN to MODIFY N_FISHES_PER_SCHOOL and N_SCHOOLS while maintaining the same N_FISHES

# 20.000 totale
for ((i = 4000; i <= 20000; i=i+4000)); do
    mpirun.actual -n 5 $EXECUTABLE_PATH_AND_NAME $i 1000 100 1    
done

for ((i = 2000; i <= 10000; i=i+2000)); do
    mpirun.actual -n 10 $EXECUTABLE_PATH_AND_NAME $i 1000 100 1    
done

for ((i = 1000; i <= 5000; i=i+1000)); do
    mpirun.actual -n 20 $EXECUTABLE_PATH_AND_NAME $i 1000 100 1    
done

for ((i = 500; i <= 2500; i=i+500)); do
    mpirun.actual -n 40 $EXECUTABLE_PATH_AND_NAME $i 1000 100 1    
done

# # 20.000 totale
# for ((i = 4000; i <= 20000; i=i+4000)); do
#     mpirun.actual -n 5 $EXECUTABLE_PATH_AND_NAME $i 1000 100 2   
# done

# for ((i = 2000; i <= 10000; i=i+2000)); do
#     mpirun.actual -n 10 $EXECUTABLE_PATH_AND_NAME $i 1000 100 2
# done

# for ((i = 1000; i <= 5000; i=i+1000)); do
#     mpirun.actual -n 20 $EXECUTABLE_PATH_AND_NAME $i 1000 100 2    
# done

# for ((i = 500; i <= 2500; i=i+500)); do
#     mpirun.actual -n 40 $EXECUTABLE_PATH_AND_NAME $i 1000 100 2   
# done

