#!/bin/bash
# max walltime 6h
#PBS -q short_cpuQ

# expected timespan for execution
#PBS -l walltime=06:00:00

# chunks (~nodes) : cores per chunk : shared memory per chunk (?)
#PBS -l select=2:ncpus=8:mem=2gb
#PBS -l place=scatter

export OMP_NUM_THREADS=8
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_SCHEDULE="dynamic, 100"

readonly C_PROGRAM_PATH=~/hpc_project/our_versions/parallel_execution/src/v3_mpi_openmp/mpi_openmp_FSO.c
readonly EXECUTABLE_PATH_AND_NAME=~/hpc_project/our_versions/parallel_execution/src/v3_mpi_openmp/mpi_openmp_FSO

# get dependencies
module load mpich-3.2

# build
mpicc $C_PROGRAM_PATH -g -Wall -fopenmp -lm -std=c99 -o $EXECUTABLE_PATH_AND_NAME 

# run
mpirun.actual -n 2 $EXECUTABLE_PATH_AND_NAME 64000 1000 200 1
