import json
import glob

#prendi il file job_ids.json e prendi gli id e le informazioni sui job
with open("results.json", "r") as f:
    job_ids = json.load(f)
    job_ids = list(job_ids.keys())

print(job_ids)

#prendi tutti i file "generated_job_*.sh.o{job_id}" in cui la 
for job_id in job_ids:
    #regex per trovare il file "generated_job_*.sh.o{job_id}" in cui la stellina sono delle info che non ci interessano
    matching_file = glob.glob(f"generated_job_*.sh.o{job_id}")
    if matching_file.__len__() >= 1:
        matching_file = matching_file[0]
    else:
        print(f"File not found for job id {job_id}")
        continue

    # trova all'interno del file una delle righe che inizia con "END" e prendi il valore successivo. 
    # Questo valore è il tempo (in secondi) che il job ha impiegato a finire.
    with open(matching_file, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("END"):
                time = line.split()[1]
                break

# Prendi questo valore e mettilo nello stesso json file di prima sotto la chiave del job_id
    with open("results.json", "r") as f:
        job_ids = json.load(f)
        job_ids[job_id]["time"] = time
    
    with open("results.json", "w") as f:
        json.dump(job_ids, f, indent=4)
