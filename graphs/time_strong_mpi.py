import json
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
        if total_fishes not in [128000, 64000, 32000, 16000]:
            continue
        nodes = values["nodes"]
        time = values["time"]
        grouped_data[total_fishes].append((nodes, time))

# Plotting
plt.figure(figsize=(8, 8))

# Ordina per numero di pesci decrescente
for total_fishes in sorted(grouped_data.keys(), reverse=True):
    points = grouped_data[total_fishes]
    points.sort(key=lambda x: x[0])
    cores, times = zip(*points)
    plt.plot(cores, times, marker='o', label=f'{total_fishes} fishes')


# plt.ylim(0, 1.3)
plt.xticks([1, 2, 4, 8, 16, 32, 64])
plt.grid(True, which='both', axis='x')
plt.xlabel('Number of Processess')
plt.ylabel('Time (s)')
plt.title(f"Time vs. Processes (places: {place_filter})")
plt.legend(title='Total Fishes')
plt.grid(True)
plt.tight_layout()
# plt.show()
plt.savefig(f"images/time_vs_nodes_mpi_strong_{place_filter}.png", dpi=100)
