import json
import matplotlib.pyplot as plt
import os

# Load JSON data from file
with open("../implementations/hybrid/results_modified.json", "r") as f:
    data = json.load(f)

# Organize data
grouped = {}

for key, entry in data.items():
    density = entry["total_fishes"] / (entry["cores"] * entry["nodes"])
    group_key = (density, entry["places"])
    if group_key not in grouped:
        grouped[group_key] = []
    grouped[group_key].append(entry)

# Show available groups to choose from
print("Available groups (density, places):")
for i, k in enumerate(grouped.keys()):
    print(f"{i}: density={k[0]:.2f}, place={k[1]}, items={len(grouped.get(k))}")

# Ask user for group selection
choice = int(input("Select a group index to plot: "))
selected_key = list(grouped.keys())[choice]
selected_data = grouped[selected_key]

# Sort by nodes
selected_data.sort(key=lambda x: (x["nodes"], x["cores"]))

# Organize by core count
from collections import defaultdict
cores_group = defaultdict(list)
for entry in selected_data:
    cores_group[entry["cores"]].append((entry["nodes"], entry["time"]))

# Plot
plt.figure(figsize=(10, 6))
for core, values in cores_group.items():
    values.sort()
    nodes, times = zip(*values)
    plt.plot(nodes, times, marker='o', label=f"{core} core(s)")

plt.xlabel("Nodes")
plt.ylabel("Time")
plt.title(f"Time vs Nodes (density={selected_key[0]:.2f}, place={selected_key[1]})")
plt.legend()
plt.grid(True)
# plt.ylim(0, 2000)  # Efficiency scale
# plt.xscale("log", base=2)  # Optional: log-scale if nodes grow exponentially
plt.xticks([1, 2, 4, 8, 16, 32, 64], labels=[1, 2, 4, 8, 16, 32, 64])  # Fixed x grid
plt.axhline(y=1.0, color='gray', linestyle='--', linewidth=1)

plt.tight_layout()
# plt.savefig(f"../images/time_vs_nodes_hybrid_weak_{selected_key[0]:.2f}_{selected_key[1]}.png", dpi=300)
plt.show()
