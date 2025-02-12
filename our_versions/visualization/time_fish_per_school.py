import matplotlib.pyplot as plt
import re

# Read log data from file
with open("../sequential_execution/src/v1.1/test.txt", "r") as file:
    log_data = file.read()

# Regular expression to extract fish count and time
pattern = re.compile(r"N_FISHES_PER_SCHOOL (\d+).*\n([\d.]+)")

# Extracting data
fish_counts = []
times = []
for match in pattern.finditer(log_data):
    fish_counts.append(int(match.group(1)))
    times.append(float(match.group(2)))

# Check if data was extracted
if not fish_counts or not times:
    print("No data extracted. Check the log format.")
else:
    # Plotting
    plt.figure(figsize=(10, 6))
    plt.plot(fish_counts, times, marker='o', linestyle='-', color='b', label="Time (ms)")
    plt.xlabel("Number of Fishes per School")
    plt.ylabel("Time (ms)")
    plt.title("Computation Time vs Number of Fishes per School")
    plt.grid(True)
    plt.legend()
    plt.show()
