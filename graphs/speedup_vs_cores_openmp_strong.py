import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("TkAgg")
from collections import defaultdict

# Load your JSON data
with open('../implementations/hybrid/results_modified.json') as f:
    data = json.load(f)


# Group data by total_fishes
grouped_data = defaultdict(list)

for run_id, values in data.items():
    if values["nodes"] == 1:
        total_fishes = values["total_fishes"]
        cores = values["cores"]
        speedup = values["speedup"]
        grouped_data[total_fishes].append((cores, speedup))

# Plotting
plt.figure(figsize=(8, 8))

for total_fishes, points in grouped_data.items():
    # Sort by number of cores for consistent lines
    points.sort(key=lambda x: x[0])
    cores, speedups = zip(*points)
    plt.plot(cores, speedups, marker='o', label=f'{total_fishes} fishes')


max_val = 75
plt.plot([0, max_val], [0, max_val], 'k--', label="linear speedup")

plt.xticks([1, 2, 4, 8, 16, 32, 64])
plt.grid(True, which='both', axis='x')
plt.xlabel('Number of Nodes')
plt.ylabel('Speedup')
plt.title("Speedup vs. Cores")
plt.legend(title='Total Fishes')
plt.grid(True)
plt.tight_layout()
# plt.show()
plt.savefig("images/speedup_vs_nodes_openmp_strong.png", dpi=100)
