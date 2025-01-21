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
#PBS -l walltime=01:30:00

# chunks (~nodes) : cores per chunk : shared memory per chunk (?)
#PBS -l select=4:ncpus=3:mem=2gb

readonly C_PROGRAM_PATH=~/hpc_project/...
readonly EXECUTABLE_PATH_AND_NAME=~/hpc_project/...

# get dependencies
module load mpich-3.2
# build
mpicc $C_PROGRAM_PATH -g -Wall -std=c99 -o $EXECUTABLE_PATH_AND_NAME
# run
mpirun.actual -n 12 $EXECUTABLE_PATH_AND_NAME
