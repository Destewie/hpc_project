import json
import numpy as np
import matplotlib
matplotlib.use("TkAgg")  # Usa il backend TkAgg per evitare problemi con Qt/xcb
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Load data from JSON file
def load_data(json_file):
    with open(json_file, 'r') as f:
        return json.load(f)

# Plot 3D scatter plot
def plot_3d_scatter(data):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    x = []
    y = []
    z = []
    colors = []

    for entry in data['values']:
        x.append(entry['number_of_fishes'])
        y.append(entry['number_of_dimensions'])
        z.append(entry['tempo_s'])
        colors.append(entry['tempo_s'])

    norm = plt.Normalize(min(colors), max(colors))
    cmap = plt.get_cmap('viridis')
    scatter_colors = cmap(norm(colors))

    sc = ax.scatter(x, y, z, c=colors, cmap='viridis', s=50)

    ax.set_xlabel('Number of Fishes')
    ax.set_ylabel('Number of Dimensions')
    ax.set_zlabel('Tempo (s)')

    plt.colorbar(sc, ax=ax, label='Tempo (s)')
    plt.show()

# Load data and plot
if __name__ == "__main__":
    data = load_data('performances.json')
    plot_3d_scatter(data)
