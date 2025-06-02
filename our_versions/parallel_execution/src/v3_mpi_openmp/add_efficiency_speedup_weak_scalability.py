import json
from collections import defaultdict

# Load your JSON data
with open("results_weak.json", "r") as f:
    data = json.load(f)

# Step 1: Organize entries by config (weak scaling: exclude total_fishes)
groups = defaultdict(list)
for key, entry in data.items():
    config = (
        entry["places"],
        entry["dimensions"],
        entry["iterations"],
        entry["update_frequency"]
    )
    groups[config].append((key, entry))

# Step 2: Process each group
for config, entries in groups.items():
    # Find the baseline with nodes == 1 and cores == 1 and smallest total_fishes
    baseline_time = None
    baseline_fishes = None
    for key, entry in entries:
        if entry["nodes"] == 1 and entry["cores"] == 1:
            if baseline_fishes is None or entry["total_fishes"] < baseline_fishes:
                baseline_time = entry["time"]
                baseline_fishes = entry["total_fishes"]

    if baseline_time is None:
        continue  # No baseline found

    # Compute speedup and efficiency for each entry
    for key, entry in entries:
        time = entry["time"]
        nodes = entry["nodes"]
        cores = entry["cores"]
        speedup = baseline_time / time
        efficiency = speedup / (nodes * cores) if nodes > 0 and cores > 0 else 0.0

        # Update entry
        entry["speedup"] = speedup
        entry["efficiency"] = efficiency

# Save modified JSON
with open("results_weak_modified.json", "w") as f:
    json.dump(data, f, indent=4)
