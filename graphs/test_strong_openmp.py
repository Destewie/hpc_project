import json
import matplotlib.pyplot as plt
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
    if values["nodes"] == 1 and values["places"] == place_filter:
        total_fishes = values["total_fishes"]
        cores = values["cores"]
        efficiency = values["efficiency"]
        grouped_data[total_fishes].append((cores, efficiency))

# Plotting
plt.figure(figsize=(10, 6))

for total_fishes, points in grouped_data.items():
    # Sort by number of cores for consistent lines
    points.sort(key=lambda x: x[0])
    cores, efficiencies = zip(*points)
    plt.plot(cores, efficiencies, marker='o', label=f'{total_fishes} fishes')

plt.xlabel('Number of Cores')
plt.ylabel('Efficiency')
plt.title(f"Efficiency vs. Cores (places: {place_filter})")
plt.legend(title='Total Fishes')
plt.grid(True)
plt.tight_layout()
plt.show()
