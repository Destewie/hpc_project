import json
import glob
import os
import re

# Percorso del file JSON
json_file = 'results_mpi_weak.json'

# Carica o crea il file JSON
if os.path.exists(json_file):
    with open(json_file, 'r') as f:
        results = json.load(f)
else:
    results = {}

# Estrai job_id già presenti
existing_ids = set(results.keys())

# Cerca file che corrispondono al pattern
all_files = glob.glob("mpi_openmp_FSO.sh.o*")

# Filtra solo quelli con ID non già presenti
new_files = []
for f in all_files:
    match = re.search(r"\.sh\.o(\d+)$", f)
    if match:
        job_id = match.group(1)
        if job_id not in existing_ids:
            new_files.append((f, job_id))

# Analizza i nuovi file
for file_path, job_id in new_files:
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Estrai i parametri da RUNNING WITH
    params_line = next((line for line in lines if line.startswith("RUNNING WITH:")), None)
    if not params_line:
        print(f"Il job {job_id} non ha una riga RUNNING WITH valida.")
        continue  # se non c'è, salta

    params_match = re.search(
        r"N-PROCESSES (\d+) - N_FISHES_PER_PROCESS (\d+) - N_THREADS (\d+) - DIMENSIONS (\d+) - MAX_ITER (\d+) - UPDATE_FREQUENCY (\d+) - PLACE (\w+)",
        params_line
    )
    if not params_match:
        print (f"Il job {job_id} non ha parametri validi.")
        continue

    nodes = int(params_match.group(1))
    fishes_per_proc = int(params_match.group(2))
    total_fishes = nodes * fishes_per_proc
    cores = int(params_match.group(3))
    dimensions = int(params_match.group(4))
    iterations = int(params_match.group(5))
    update_freq = int(params_match.group(6))
    place = params_match.group(7)

    # Cerca riga END e prendine un valore con parte decimale
    time = None
    for line in lines:
        if line.startswith("END:"):
            parts = line.strip().split(":")
            if len(parts) > 1:
                try:
                    time = float(parts[1])
                    break
                except ValueError:
                    print(f"Il job {job_id} ha un tempo non valido: {parts[1]}")
                    time = None
                    break


    if time is None:
        print(f"Il job {job_id} non ha un tempo di esecuzione valido.")
        continue

    # Salva i risultati
    results[job_id] = {
        "total_fishes": total_fishes,
        "dimensions": dimensions,
        "iterations": iterations,
        "update_frequency": update_freq,
        "nodes": nodes,
        "cores": cores,
        "places": place,  # fissa come nel tuo esempio
        "time": time
    }

# Scrive su disco il file aggiornato
with open(json_file, 'w') as f:
    json.dump(results, f, indent=4)

print(f"{len(new_files)} file analizzati e salvati.")
