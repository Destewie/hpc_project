import matplotlib.pyplot as plt
import re

log_data = "../sequential_execution/v1.1/sequential_FSO.sh.o2710978"  # Replace this with your actual log data

# Regular expression to extract data
pattern = re.compile(r"N_FISHES_PER_SCHOOL (\d+).+?\n(\d+\.\d+)")

# Extracting data
fish_counts = []
times = []
for match in pattern.finditer(log_data):
    fish_counts.append(int(match.group(1)))
    times.append(float(match.group(2)))

# Plotting
title = "Computation Time vs Number of Fishes per School"
plt.figure(figsize=(10, 6))
plt.plot(fish_counts, times, marker='o', linestyle='-', color='b', label="Time (ms)")
plt.xlabel("Number of Fishes per School")
plt.ylabel("Time (ms)")
plt.title(title)
plt.grid(True)
plt.legend()
plt.show()