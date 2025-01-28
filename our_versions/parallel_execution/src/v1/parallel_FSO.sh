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
#PBS -l walltime=02:00:00

# chunks (~nodes) : cores per chunk : shared memory per chunk (?)
#PBS -l select=2:ncpus=5:mem=2gb

readonly C_PROGRAM_PATH=~/hpc_project/our_versions/parallel_execution/src/v1/parallell_FSO.c
readonly EXECUTABLE_PATH_AND_NAME=~/hpc_project/our_versions/parallel_execution/src/v1/parallell_FSO

# get dependencies
module load mpich-3.2
# build
mpicc $C_PROGRAM_PATH -g -Wall -fopenmp -lm -std=c99 -o $EXECUTABLE_PATH_AND_NAME 
# run
mpirun.actual -n 10 $EXECUTABLE_PATH_AND_NAME
