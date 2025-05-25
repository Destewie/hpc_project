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
#PBS -l select=4:ncpus=4:mem=2gb
#PBS -l place=scatter

export OMP_NUM_THREADS=4 
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_SCHEDULE="dynamic, 1"


readonly C_PROGRAM_PATH=~/hpc_project/our_versions/parallel_execution/src/v3_mpi_openmp/mpi_openmp_FSO.c
readonly EXECUTABLE_PATH_AND_NAME=~/hpc_project/our_versions/parallel_execution/src/v3_mpi_openmp/mpi_openmp_FSO

# get dependencies
module load mpich-3.2

# build
mpicc $C_PROGRAM_PATH -g -Wall -fopenmp -lm -std=c99 -o $EXECUTABLE_PATH_AND_NAME 

# run
# <"Usage: N_FISHES_PER_PROCESS DIMENSIONS MAX_ITER UPDATE_FREQUENCY">
mpirun.actual -n 4 $EXECUTABLE_PATH_AND_NAME 16000 1000 200 1 scatter


