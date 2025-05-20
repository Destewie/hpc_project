import os
import subprocess

# VALID_SELECT = [1, 2, 4, 8, 16, 32, 64]
# VALID_NCPUS = [1, 2, 4, 8, 16, 32]

VALID_SELECT = [2]
VALID_NCPUS = [4,8]
VALID_PLACE = ['pack', 'scatter']
TOTAL_FISHES = 128000

PBS_TEMPLATE = """#!/bin/bash
# max walltime 6h
#PBS -q short_cpuQ

# expected timespan for execution
#PBS -l walltime=06:00:00

# chunks (~nodes) : cores per chunk : shared memory per chunk (?)
#PBS -l select={select}:ncpus={ncpus}:mem=2gb
#PBS -l place={place}

export OMP_NUM_THREADS={ncpus}
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
mpirun.actual -n {select} $EXECUTABLE_PATH_AND_NAME {fishes_per_process} 1000 200 1
"""

def generate_pbs_script(select, ncpus, place, output_path="generated_job.sh"):
    fishes_per_process = TOTAL_FISHES/select
    pbs_script = PBS_TEMPLATE.format(select=select, ncpus=ncpus, place=place, fishes_per_process=fishes_per_process)

    with open(output_path, "w") as f:
        f.write(pbs_script)
    
    print("PBS script generated at ", output_path)

if __name__ == "__main__":
    for node in VALID_SELECT:
        for core in VALID_NCPUS:
            for place in VALID_PLACE:   
                generate_pbs_script(node, core, place)
                print("running with ", node, "processes, ", core, "cores, place=", place)
                subprocess.call(["qsub", "generated_job.sh"])
                print("\n")