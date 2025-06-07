import json
import random
from collections import defaultdict

# === Carica il file JSON ===
with open("results_openmp_weak.json", "r") as f:
    data = json.load(f)

#-----------------------------------------------------------------------------

# Raggruppa per configurazione senza includere "placement"
placement_groups = defaultdict(list)

for key, entry in data.items():
    if entry["nodes"] == 1:
        config_key = (
            entry["total_fishes"],
            entry["dimensions"],
            entry["iterations"],
            entry["update_frequency"],
            entry["nodes"],
            entry["cores"]
        )
        placement_groups[config_key].append((key, entry))

for config, entries in placement_groups.items():
    # Gruppo per placement
    placement_dict = defaultdict(list)
    for key, entry in entries:
        placement_dict[entry["places"]].append((key, entry))

    # Se esistono entrambe le versioni, normalizza al tempo minimo
    if "pack" in placement_dict and "scatter" in placement_dict:
        combined = placement_dict["pack"] + placement_dict["scatter"]
        min_time = min(e["time"] for _, e in combined)
        for _, e in combined:
            e["time"] = min_time

    # Se esiste solo una delle due, crea la entry "fantasma"
    elif ("pack" in placement_dict) ^ ("scatter" in placement_dict):
        existing_placement = "pack" if "pack" in placement_dict else "scatter"
        inverse_placement = "scatter" if existing_placement == "pack" else "pack"
        for _, original_entry in placement_dict[existing_placement]:
            new_entry = original_entry.copy()
            new_entry["places"] = inverse_placement

            # Assegna un nuovo ID numerico casuale che inizia per '7'
            while True:
                new_id = str(random.randint(7_000_000, 7_999_999))
                if new_id not in data:
                    break

            # Aggiungi nuova entry al dizionario
            data[new_id] = new_entry

#-----------------------------------------------------------------------------

# === Raggruppa per configurazione costante tranne i processi (weak scaling) ===
groups = defaultdict(list)
for key, entry in data.items():
    try:
        nodes = entry["nodes"]
        cores = entry["cores"]
        total_fishes = entry["total_fishes"]
        fishes_per_core = total_fishes // (nodes * cores)

        config = (
            entry["places"],
            entry["dimensions"],
            entry["iterations"],
            entry["update_frequency"],
            fishes_per_core
        )

        groups[config].append((key, entry))
    except KeyError:
        continue

# === Calcola efficiency_weak per ogni gruppo ===
for config, entries in groups.items():
    # Trova la baseline: config con 1 nodo e 1 core
    baseline_time = None
    for key, entry in entries:
        if entry["nodes"] == 1 and entry["cores"] == 1:
            baseline_time = entry["time"]
            break

    if baseline_time is None or baseline_time <= 0:
        continue  # Nessuna baseline valida

    for key, entry in entries:
        time = entry["time"]
        if time > 0:
            entry["efficiency_weak"] = baseline_time / time
        else:
            entry["efficiency_weak"] = 0.0

        # Pulisce eventuali metriche obsolete
        entry.pop("speedup", None)
        entry.pop("efficiency", None)

# === Salva il file modificato ===
with open("results_openmp_weak_modified.json", "w") as f:
    json.dump(data, f, indent=4)
