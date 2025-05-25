import json
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
import numpy as np
import matplotlib
matplotlib.use("TkAgg")


# === Leggi il file JSON ===
with open("results.json", "r") as f:
    data = json.load(f)

# === Raccogli tutti i valori unici di total_fishes ===
fishes_options = sorted({entry["total_fishes"] for entry in data.values()})
print("Valori disponibili per 'total_fishes':")
for i, val in enumerate(fishes_options):
    print(f"[{i}] {val}")

# === Seleziona uno dei valori ===
index = int(input("Seleziona l'indice del valore di 'total_fishes' da visualizzare: "))
selected_fishes = fishes_options[index]

# === Filtra i dati ===
filtered = [entry for entry in data.values() if entry["total_fishes"] == selected_fishes]

# === Estrai valori per il grafico ===
nodes = [entry["nodes"] for entry in filtered]
cores = [entry["cores"] for entry in filtered]
times = [entry["time"] for entry in filtered]

# === Normalizza i tempi per la colormap ===
norm = plt.Normalize(min(times), max(times))
colors = cm.get_cmap('RdYlGn_r')(norm(times))  # scala: rosso (alto) → verde (basso)

# === Crea il grafico ===
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

sc = ax.scatter(nodes, cores, times, c=colors, s=60)

ax.set_xlabel("Nodes")
ax.set_ylabel("Cores")
ax.set_zlabel("Time (s)")
ax.set_title(f"Performance per total_fishes = {selected_fishes}")

# Aggiungi barra colore
mappable = cm.ScalarMappable(norm=norm, cmap='RdYlGn_r')
mappable.set_array(times)
fig.colorbar(mappable, label="Time (s)")

plt.show()
