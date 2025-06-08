import json
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib
matplotlib.use("TkAgg")
import numpy as np
import argparse

# Costante modificabile
TARGET_TOTAL_CORES = 16

# CLI args
parser = argparse.ArgumentParser(description="Plot execution time for different (nodes, cores) combinations.")
parser.add_argument("--places", choices=["scatter", "pack"], required=True, help="Placement strategy to filter by.")
args = parser.parse_args()
selected_places = args.places

# Carica i dati
with open("../implementations/hybrid/results_modified.json") as f:
    data = json.load(f)

# Genera combinazioni (nodes, cores) tali che nodes * cores == TARGET_TOTAL_CORES
valid_combinations = [
    (n, TARGET_TOTAL_CORES // n)
    for n in sorted([i for i in range(1, TARGET_TOTAL_CORES + 1) if TARGET_TOTAL_CORES % i == 0], reverse=True)
]

# Filtro delle entry
combo_to_time = {}
combo_to_eff = {}
for (n, c) in valid_combinations:
    matching_entries = []
    for entry in data.values():
        if (
            entry.get("total_fishes") == 64000 and
            entry.get("dimensions") == 1000 and
            entry.get("iterations") == 200 and
            entry.get("update_frequency") == 1 and
            entry.get("nodes") == n and
            entry.get("cores") == c
        ):
            if n == 1:
                matching_entries.append(entry)  # Ignora 'places' se n == 1
            elif entry.get("places") == selected_places:
                matching_entries.append(entry)

    if matching_entries:
        best_time = min(e["time"] for e in matching_entries)
        best_eff = max(e["efficiency"] for e in matching_entries)
        combo_to_time[(n, c)] = best_time
        combo_to_eff[(n, c)] = best_eff

# Prepara i dati per il grafico
labels = [f"({n},{c})" for (n, c) in valid_combinations if (n, c) in combo_to_time]
times = [combo_to_time[(n, c)] for (n, c) in valid_combinations if (n, c) in combo_to_time]
efficiencies = [combo_to_eff[(n, c)] for (n, c) in valid_combinations if (n, c) in combo_to_eff]
for i in range(len(times)):
    #limit times to 2 digits

    print (labels[i], "& ", round(times[i],3), " & ", round(efficiencies[i],3), "\\\\")

# Color mapping (senza colorbar)
norm = plt.Normalize(min(times), max(times))
colors = [plt.cm.RdYlGn_r(norm(t)) for t in times]

# Plot
fig, ax = plt.subplots(figsize=(10, 6))
bars = ax.bar(labels, times, color=colors)
ax.set_xlabel("(Nodes, Cores)")
ax.set_ylabel("Time (s)")
ax.set_title(f"Execution Time for (Nodes x Cores = {TARGET_TOTAL_CORES}), placement = {selected_places}")
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, rotation=45)

plt.ylim(0, 1400)
plt.tight_layout()
plt.show()
# plt.savefig(f"images/compared_times_{selected_places}_{TARGET_TOTAL_CORES}.png", dpi=100)
