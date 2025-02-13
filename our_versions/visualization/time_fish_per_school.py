import matplotlib.pyplot as plt
import re

# Read log data from file
with open("../sequential_execution/src/v1.1/same_nfish2.txt", "r") as file:
    log_data = file.read()

# Regular expression to extract data
pattern = re.compile(r"N-SCHOOLS (\d+) - N_FISHES_PER_SCHOOL (\d+).+?\n(\d+\.\d+)")

# Extracting data
data = {}
for match in pattern.finditer(log_data):
    n_schools = int(match.group(1))
    fish_count = int(match.group(2))
    time = float(match.group(3))/1000.0
    
    if n_schools not in data:
        data[n_schools] = {"fish_counts": [], "times": []}
    data[n_schools]["fish_counts"].append(fish_count)
    data[n_schools]["times"].append(time)

# Plotting
title = "Computation Time vs Number of Fishes per School"
plt.figure(figsize=(10, 6))
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']  # Define a set of colors

for i, (n_schools, values) in enumerate(sorted(data.items())):
    plt.plot(values["fish_counts"], values["times"], marker='o', linestyle='-', color=colors[i % len(colors)], label=f"N-SCHOOLS {n_schools}")

plt.xlabel("Number of Fishes per School")
plt.ylabel("Time (s)")
plt.title(title)
plt.grid(True)
plt.legend()
plt.show()
