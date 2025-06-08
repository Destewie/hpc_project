import json
from collections import defaultdict

# --- Configuration ---
target_place = "pack"  # or "pack"
target_cores = 1
fish_order = [128000, 64000, 32000, 16000]
metrics = ["time", "speedup", "efficiency"]  # will print one table per metric

# --- Load JSON ---
with open("../implementations/hybrid/results_modified.json", "r") as f:
    data = json.load(f)

# --- Collect data: {metric -> {nodes -> {total_fishes -> value}}} ---
tables = {metric: defaultdict(dict) for metric in metrics}

for entry in data.values():
    if entry["places"] != target_place:
        continue
    if entry["cores"] != target_cores:
        continue
    fishes = entry["total_fishes"]
    if fishes not in fish_order:
        continue
    nodes = entry["nodes"]
    for metric in metrics:
        tables[metric][nodes][fishes] = entry[metric]

# --- Format tables ---
for metric in metrics:
    print(f"\n### {metric.upper()} TABLE ###")
    header = ["Processes"] + [str(fish) for fish in fish_order]
    print(" & ".join(header))

    for nodes in sorted(tables[metric].keys()):
        row = [str(nodes)]
        for fish in fish_order:
            value = tables[metric][nodes].get(fish, "N/A")
            row.append(f"{value:.4f}" if isinstance(value, float) else str(value))
        print(" & ".join(row))
