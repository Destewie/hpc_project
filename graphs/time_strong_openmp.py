import json
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
        time = values["time"]
        grouped_data[total_fishes].append((cores, time))

# Plotting
plt.figure(figsize=(8, 8))

# Ordina per numero di pesci decrescente
for total_fishes in sorted(grouped_data.keys(), reverse=True):
    if total_fishes in [128000, 64000, 32000, 16000]:
        points = grouped_data[total_fishes]
        points.sort(key=lambda x: x[0])
        cores, times = zip(*points)
        plt.plot(cores, times, marker='o', label=f'{total_fishes} fishes')


# plt.ylim(0, 1.3)
plt.xticks([1, 2, 4, 8, 16, 32, 64])
plt.grid(True, which='both', axis='x')
plt.xlabel('Number of Threads (single process)')
plt.ylabel('Time (s)')
plt.title("Time vs. Threads (single process)")
plt.legend(title='Total Fishes')
plt.grid(True)
plt.tight_layout()
# plt.show()
plt.savefig(f"images/time_vs_cores_openmp_strong.png", dpi=100)
