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
#PBS -l select=1:ncpus=4:mem=2gb
export OMP_NUM_THREADS=4
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_SCHEDULE="dynamic, 100"


readonly C_PROGRAM_PATH=~/hpc_project/our_versions/parallel_execution/src/v2_multithreading_only/multithreading_FSO.c
readonly EXECUTABLE_PATH_AND_NAME=~/hpc_project/our_versions/parallel_execution/src/v2_multithreading_only/multithreading_FSO

# get dependencies
module load mpich-3.2

# build
mpicc $C_PROGRAM_PATH -g -Wall -fopenmp -lm -std=c99 -o $EXECUTABLE_PATH_AND_NAME 

# run
# <"Usage: N_SCHOOLS N_FISHES_PER_SCHOOL DIMENSIONS MAX_ITER UPDATE_FREQUENCY">
mpirun.actual -n 1 $EXECUTABLE_PATH_AND_NAME 1 5000 1000 100 1

# PLAN to MODIFY N_FISHES_PER_SCHOOL and N_SCHOOLS while maintaining the same N_FISHES
# -n x -> x is the number of processes = number of schools
# first param   = number of fishes per school
# second param  = number of dimensions
# third param   = number of iterations
# fourth param  = update frequency 

# # 20.000 pesci totali - update freq = 1
# for ((i = 4000; i <= 20000; i=i+4000)); do
#     mpirun.actual -n 5 $EXECUTABLE_PATH_AND_NAME $i 1000 100 1    
# done

# for ((i = 2000; i <= 10000; i=i+2000)); do
#     mpirun.actual -n 10 $EXECUTABLE_PATH_AND_NAME $i 1000 100 1    
# done

# for ((i = 1000; i <= 5000; i=i+1000)); do
#     mpirun.actual -n 20 $EXECUTABLE_PATH_AND_NAME $i 1000 100 1    
# done

# for ((i = 500; i <= 2500; i=i+500)); do
#     mpirun.actual -n 40 $EXECUTABLE_PATH_AND_NAME $i 1000 100 1    
# done

# # 20.000 pesci totali - update freq = 2
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

