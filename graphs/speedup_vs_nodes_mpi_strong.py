import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("TkAgg")
from collections import defaultdict

# Load your JSON data
with open('../implementations/hybrid/results_modified.json') as f:
    data = json.load(f)

# Ask user for which "places" type to show
choice = input("Enter 0 for 'pack', 1 for 'scatter': ")
place_filter = 'pack' if choice == '0' else 'scatter'

# Group data by total_fishes
grouped_data = defaultdict(list)

for run_id, values in data.items():
    if values["cores"] == 1 and values["places"] == place_filter:
        total_fishes = values["total_fishes"]
        nodes = values["nodes"]
        speedup = values["speedup"]
        grouped_data[total_fishes].append((nodes, speedup))

# Plotting
plt.figure(figsize=(8, 8))

for total_fishes, points in grouped_data.items():
    if total_fishes in [128000, 64000, 32000, 16000]:
        # Sort by number of cores for consistent lines
        points.sort(key=lambda x: x[0])
        cores, efficiencies = zip(*points)
        plt.plot(cores, efficiencies, marker='o', label=f'{total_fishes} fishes')


max_val = max(64,75)
plt.plot([0, max_val], [0, max_val], 'k--', label="linear speedup")

plt.ylim(0, 75)
plt.xticks([1, 2, 4, 8, 16, 32, 64])
plt.grid(True, which='both', axis='x')
plt.xlabel('Number of Processes')
plt.ylabel('Speedup')
plt.title(f"Speedup vs. Processes (places: {place_filter})")
plt.legend(title='Total Fishes')
plt.grid(True)
plt.tight_layout()
# plt.show()
plt.savefig(f"images/speedup_vs_nodes_mpi_strong_{place_filter}.png", dpi=100)
