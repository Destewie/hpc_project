import json
from collections import defaultdict

# Load your JSON data
with open("results.json", "r") as f:
    data = json.load(f)

# Step 1: Organize entries by problem configuration
groups = defaultdict(list)
for key, entry in data.items():
    config = (
        entry["total_fishes"],
        entry["dimensions"],
        entry["iterations"],
        entry["update_frequency"]
    )
    groups[config].append((key, entry))

# Step 2: Process each group
for config, entries in groups.items():
    # Find the baseline entry: nodes == 1 and cores == 1
    baseline_time = None
    for key, entry in entries:
        if entry["nodes"] == 1 and entry["cores"] == 1:
            baseline_time = entry["time"]
            break

    # If no baseline, skip this group
    if baseline_time is None:
        continue

    # Step 3: Compute speedup and efficiency
    for key, entry in entries:
        time = entry["time"]
        nodes = entry["nodes"]
        speedup = baseline_time / time
        efficiency = speedup / nodes if nodes > 0 else 0.0

        # Update entry
        entry["speedup"] = speedup
        entry["efficiency"] = efficiency

# Save modified JSON
with open("results_modified.json", "w") as f:
    json.dump(data, f, indent=4)
