import json
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from collections import defaultdict
import sys

# === Carica il file JSON ===
with open('../implementations/openmp_only/results_openmp_weak_modified.json', 'r') as f:
    data = json.load(f)

# === Raggruppa per fishes_per_core e config costante ===
groups = defaultdict(list)

for key, entry in data.items():
    try:
        cores = entry["cores"]
        nodes = entry["nodes"]
        total_fishes = entry["total_fishes"]
        fishes_per_core = total_fishes // (cores * nodes)

        group_key = (
            fishes_per_core,
            entry["dimensions"],
            entry["iterations"],
            entry["places"],
            entry["update_frequency"]
        )
        groups[group_key].append((key, entry))
    except KeyError:
        continue

# === Filtra gruppi validi ===
valid_groups = {}

for group_key, entries in groups.items():
    if len(entries) < 2:  # almeno 2 entry per avere un grafico utile
        continue

    iterations_set = {e["iterations"] for _, e in entries}
    dimensions_set = {e["dimensions"] for _, e in entries}
    places_set = {e["places"] for _, e in entries}

    if len(iterations_set) == 1 and len(dimensions_set) == 1 and len(places_set) == 1:
        valid_groups[group_key] = entries

if not valid_groups:
    print("Nessun gruppo valido trovato con almeno 2 elementi.")
    sys.exit()

# === Mostra i gruppi e chiedi input CLI ===
print("Gruppi validi:")
group_keys = list(valid_groups.keys())
for i, group_key in enumerate(group_keys):
    fishes_per_core, dimensions, iterations, places, freq = group_key
    print(f"[{i}] Fishes/core: {fishes_per_core}, Dim: {dimensions}, Iter: {iterations}, Places: {places}, Freq: {freq}, Entries: {len(valid_groups[group_key])}")

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
    efficiency_weak = entry.get("efficiency_weak", None)
    if efficiency_weak is not None:
        x.append(total_cores)
        y.append(efficiency_weak)

# === Plot ===
plt.figure(figsize=(8, 8))
plt.axhline(y=1.0, color='gray', linestyle='--', linewidth=1, label='Ideal Efficiency')
plt.ylim(0, 1.3)
plt.xticks([1, 2, 4, 8, 16, 32, 64])
plt.grid(True, which='both', axis='x')
plt.plot(x, y, marker='o', linestyle='-', color='darkorange')
plt.title(f"Weak Scalability Efficiency vs Total Cores")
plt.xlabel("Total Cores (nodes × cores)")
plt.ylabel("Efficiency (Weak)")
plt.xticks(x)
plt.grid(True)
plt.tight_layout()
# plt.show()
plt.savefig(f"images/efficiency_vs_nodes_openmp_weak.png", dpi=100)