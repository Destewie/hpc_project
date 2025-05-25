import json
import matplotlib
matplotlib.use("TkAgg")  # Forza il backend TkAgg per compatibilità
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# === Caricamento dati JSON ===
with open("results.json", "r") as f:
    data = json.load(f)

# === Estrai opzioni uniche per total_fishes ===
fishes_options = sorted({entry["total_fishes"] for entry in data.values()})
print("Valori disponibili per 'total_fishes':")
for i, val in enumerate(fishes_options):
    print(f"[{i}] {val}")

# === Input utente ===
index = int(input("Seleziona l'indice del valore di 'total_fishes' da visualizzare: "))
selected_fishes = fishes_options[index]

# === Filtro dati ===
filtered = [entry for entry in data.values() if entry["total_fishes"] == selected_fishes]

nodes = [entry["nodes"] for entry in filtered]
cores = [entry["cores"] for entry in filtered]
times = [entry["time"] for entry in filtered]

# === Normalizzazione colori ===
norm = plt.Normalize(min(times), max(times))

# === Seleziona colormap in modo robusto ===
try:
    cmap = matplotlib.colormaps['RdYlGn_r']  # Matplotlib >= 3.8
except AttributeError:
    from matplotlib import cm
    cmap = cm.get_cmap('RdYlGn_r')  # fallback per versioni < 3.8

colors = cmap(norm(times))

# === Plot 3D ===
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
sc = ax.scatter(nodes, cores, times, c=colors, s=60)

ax.set_xlabel("Nodes")
ax.set_ylabel("Cores")
ax.set_zlabel("Time (s)")
ax.set_title(f"Performance per total_fishes = {selected_fishes}")

# === Aggiunta colorbar correttamente ===
from matplotlib.cm import ScalarMappable
mappable = ScalarMappable(norm=norm, cmap=cmap)
mappable.set_array(times)
fig.colorbar(mappable, ax=ax, label="Time (s)")  # associato esplicitamente all'Axes

plt.show()
