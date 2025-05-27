import json
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from collections import defaultdict
import sys

# === Carica il file JSON ===
with open('results_modified.json', 'r') as f:
    data = json.load(f)

# === Raggruppa per (cores, fishes_per_core) ===
groups = defaultdict(list)

for key, entry in data.items():
    try:
        cores = entry["cores"]
        nodes = entry["nodes"]
        total_fishes = entry["total_fishes"]
        fishes_per_core = total_fishes // (cores * nodes)

        group_key = (cores, fishes_per_core)
        groups[group_key].append((key, entry))
    except KeyError:
        continue

# === Filtra gruppi validi ===
valid_groups = {}

for group_key, entries in groups.items():
    if len(entries) < 3:
        continue

    iterations_set = {e["iterations"] for _, e in entries}
    dimensions_set = {e["dimensions"] for _, e in entries}
    places_set = {e["places"] for _, e in entries}

    if len(iterations_set) == 1 and len(dimensions_set) == 1 and len(places_set) == 1:
        valid_groups[group_key] = entries

if not valid_groups:
    print("Nessun gruppo valido trovato con almeno 3 elementi.")
    sys.exit()

# === Mostra i gruppi e chiedi input CLI ===
print("Gruppi validi:")
group_keys = list(valid_groups.keys())
for i, (cores, fishes_per_core) in enumerate(group_keys):
    print(f"[{i}] Cores per proc: {cores}, Fishes/core: {fishes_per_core}, Entries: {len(valid_groups[(cores, fishes_per_core)])}")

try:
    selected_index = int(input("Seleziona il numero del gruppo da analizzare: "))
    if selected_index < 0 or selected_index >= len(group_keys):
        raise ValueError
except ValueError:
    print("Input non valido.")
    sys.exit()

# === Estrai gruppo selezionato ===
selected_key = group_keys[selected_index]
selected_group = valid_groups[selected_key]

# === Prepara dati per il grafico ===
x = []
y = []

for _, entry in sorted(selected_group, key=lambda e: e[1]["nodes"] * e[1]["cores"]):
    total_cores = entry["nodes"] * entry["cores"]
    x.append(total_cores)
    y.append(entry["efficiency"])

# === Plot ===
plt.figure(figsize=(10, 6))
plt.plot(x, y, marker='o', linestyle='-', color='darkorange')
plt.title(f"Efficiency vs Total Cores (Group {selected_index})")
plt.xlabel("Total Cores (nodes × cores)")
plt.ylabel("Efficiency")
plt.xticks(x)
plt.grid(True)
plt.tight_layout()
plt.show()
