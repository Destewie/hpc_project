import subprocess
import json

# se proviamo a lanciare tutte le combinazioni tra processi, thread e place, produciamo troppi job per il cluster
# teniamo fissi i processi e lanciamo tutte le altre combinazioni
# VALID_SELECT = [1, 2, 4, 8, 16, 32, 64]
VALID_SELECT = [4]
# VALID_NCPUS = [1, 2, 4, 8, 16, 32]
VALID_NCPUS = [1, 2]
# VALID_PLACE = ['pack', 'scatter']
VALID_PLACE = ['scatter']
TOTAL_FISHES = 16000

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
export OMP_SCHEDULE="dynamic, 1"

readonly C_PROGRAM_PATH=~/hpc_project/our_versions/parallel_execution/src/v3_mpi_openmp/mpi_openmp_FSO.c
readonly EXECUTABLE_PATH_AND_NAME=~/hpc_project/our_versions/parallel_execution/src/v3_mpi_openmp/mpi_openmp_FSO_{select}_{ncpus}_{place}.out

# get dependencies
module load mpich-3.2

# build
mpicc $C_PROGRAM_PATH -g -Wall -fopenmp -lm -std=c99 -o $EXECUTABLE_PATH_AND_NAME 

# run
mpirun.actual -n {select} $EXECUTABLE_PATH_AND_NAME {fishes_per_process} 1000 200 1 {place}
"""

def generate_pbs_script(select, ncpus, place, output_path="generated_job.sh"):
    fishes_per_process = TOTAL_FISHES/select
    pbs_script = PBS_TEMPLATE.format(select=select, ncpus=ncpus, place=place, fishes_per_process=fishes_per_process)

    output_path = f"generated_job_{select}_{ncpus}_{place}.sh" 
    with open(output_path, "w") as f:
        f.write(pbs_script)

    return output_path

    
if __name__ == "__main__":
    ids = []
    nodes = []
    cores = []
    places = []

    for node in VALID_SELECT:
        for core in VALID_NCPUS:
            for place in VALID_PLACE:   
                job_file_string = generate_pbs_script(node, core, place)

                print("\n")
                print(f"running with {node} processes, {core} cores, place={place}")
                result = subprocess.run(
                    ["qsub", job_file_string],
                    capture_output=True,
                    text=True  
                )

                # remove everything from result.stdout after the first .
                # from 3120330.hpc-head-n1.unitn.it to 3120330
                job_id = result.stdout.split(".")[0]
                print("JOB_ID:", job_id, "full:", result.stdout)
                print("STDERR:", result.stderr)
                print("Return code:", result.returncode)
                
                if result.stderr == "" and result.returncode == 0:
                    ids.append(job_id)
                    nodes.append(node)
                    cores.append(core)
                    places.append(place)

    # create a json file with the job ids as keys and the parameters with the same index as its values
    with open("job_ids.json", "w") as f:
        json.dump({ids[i]: {"nodes": nodes[i], "cores": cores[i], "places": places[i]} for i in range(len(ids))}, f, indent=4)


                
