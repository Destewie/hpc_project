import json
from collections import defaultdict

# === Carica il file JSON ===
with open("results_openmp_weak.json", "r") as f:
    data = json.load(f)

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
