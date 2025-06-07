import json
import matplotlib
matplotlib.use("TkAgg")  # Usa TkAgg per evitare problemi con Qt
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

# === Leggi il file JSON ===
with open("../implementations/hybrid/results_modified.json", "r") as f:
    data = json.load(f)

# === Trova i valori unici di total_fishes ===
fishes_options = sorted({entry["total_fishes"] for entry in data.values()})
print("Valori disponibili per 'total_fishes':")
for i, val in enumerate(fishes_options):
    print(f"[{i}] {val}")

# === Scelta utente ===
index = int(input("Seleziona l'indice del valore di 'total_fishes' da visualizzare: "))
selected_fishes = fishes_options[index]

# === Trova i valori unici per 'places' ===
places_options = sorted({entry["places"] for entry in data.values() if entry["total_fishes"] == selected_fishes})
print("\nValori disponibili per 'places':")
for i, val in enumerate(places_options):
    print(f"[{i}] {val}")

place_index = int(input("Seleziona l'indice del valore di 'places' da visualizzare: "))
selected_place = places_options[place_index]

# === Filtra i dati ===
filtered = [
    entry for entry in data.values()
    if entry["total_fishes"] == selected_fishes and entry["places"] == selected_place
]

# === Raggruppa per numero di nodi ===
grouped_by_nodes = defaultdict(list)
for entry in filtered:
    nodes = entry["nodes"]
    grouped_by_nodes[nodes].append((entry["cores"], entry["time"]))

# === Crea grafico 2D ===
plt.figure(figsize=(8, 8))
for nodes, values in sorted(grouped_by_nodes.items()):
    if nodes != 64:
        values.sort()
        x = [cores for cores, _ in values]
        y = [time for _, time in values]
        plt.plot(x, y, marker='o', label=f"{nodes} nodo{'i' if nodes > 1 else ''}")

plt.xticks([1, 2, 4, 8, 16, 32, 64])
plt.xlim(0, 67)
plt.ylim(0, 3700)  # Imposta un limite superiore per il tempo
plt.xlabel("Numero di core")
plt.ylabel("Tempo (s)")
plt.title(f"Tempo vs Core (fishes={selected_fishes}, places={selected_place})")
plt.legend(title="Numero di nodi")
plt.grid(True)
plt.tight_layout()
# plt.show()
plt.savefig(f"images/time_vs_cores_hybrid_strong_{selected_place}_{selected_fishes}.png", dpi=100)



