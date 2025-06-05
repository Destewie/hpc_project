import json
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib.cm import ScalarMappable
from scipy.interpolate import griddata

# === Caricamento dati JSON ===
with open("results.json", "r") as f:
    data = json.load(f)

# === Estrai opzioni uniche per total_fishes ===
fishes_options = sorted({entry["total_fishes"] for entry in data.values()})
print("Valori disponibili per 'total_fishes':")
for i, val in enumerate(fishes_options):
    print(f"[{i}] {val}")
index = int(input("Seleziona l'indice di 'total_fishes': "))
selected_fishes = fishes_options[index]

# === Estrai opzioni places disponibili per quel valore ===
places_options = sorted({entry["places"] for entry in data.values() if entry["total_fishes"] == selected_fishes})
print("\nValori disponibili per 'places':")
for i, val in enumerate(places_options):
    print(f"[{i}] {val}")
place_index = int(input("Seleziona l'indice di 'places': "))
selected_place = places_options[place_index]

# === Filtra dati ===
filtered = [
    entry for entry in data.values()
    if entry["total_fishes"] == selected_fishes and entry["places"] == selected_place
]

nodes = np.array([entry["nodes"] for entry in filtered])
cores = np.array([entry["cores"] for entry in filtered])
times = np.array([entry["time"] for entry in filtered])

# === Normalizzazione colori ===
norm = plt.Normalize(min(times), max(times))
try:
    cmap = matplotlib.colormaps['RdYlGn_r']
except AttributeError:
    from matplotlib import cm
    cmap = cm.get_cmap('RdYlGn_r')
colors = cmap(norm(times))

# === Preparazione figura 3D ===
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# === Superficie interpolata (piano) ===
grid_nodes, grid_cores = np.meshgrid(
    np.linspace(min(nodes), max(nodes), 30),
    np.linspace(min(cores), max(cores), 30)
)
grid_times = griddata((nodes, cores), times, (grid_nodes, grid_cores), method='cubic')

# === Visualizza superficie interpolata ===
surf = ax.plot_surface(
    grid_nodes, grid_cores, grid_times,
    cmap=cmap, norm=norm, alpha=0.7, edgecolor='k', linewidth=0.2
)

# === Scatter dei punti originali ===
sc = ax.scatter(nodes, cores, times, c=colors, s=50)

# === Label e titolo ===
ax.set_xlabel("Nodes")
ax.set_ylabel("Cores")
ax.set_zlabel("Time (s)")
ax.set_title(f"Performance per total_fishes = {selected_fishes}, places = {selected_place}")

# === Colorbar ===
mappable = ScalarMappable(norm=norm, cmap=cmap)
mappable.set_array(times)
fig.colorbar(mappable, ax=ax, label="Time (s)")

plt.tight_layout()
plt.show()
