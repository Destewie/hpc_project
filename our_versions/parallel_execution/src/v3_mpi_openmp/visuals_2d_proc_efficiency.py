import json
import matplotlib
matplotlib.use("TkAgg")  # Usa TkAgg per evitare problemi con Qt
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

# === Leggi il file JSON ===
with open("results_modified.json", "r") as f:
    data = json.load(f)

# === Trova i valori unici di total_fishes ===
fishes_options = sorted({entry["total_fishes"] for entry in data.values()})
print("Valori disponibili per 'total_fishes':")
for i, val in enumerate(fishes_options):
    print(f"[{i}] {val}")

# === Scelta utente ===
index = int(input("Seleziona l'indice del valore di 'total_fishes' da visualizzare: "))
selected_fishes = fishes_options[index]

# === Filtra i dati ===
filtered = [entry for entry in data.values() if entry["total_fishes"] == selected_fishes]

# === Raggruppa per numero di nodi ===
grouped_by_nodes = defaultdict(list)
for entry in filtered:
    nodes = entry["cores"]
    grouped_by_nodes[nodes].append((entry["nodes"], entry["speedup"]/ entry["cores"]))

# === Crea grafico 2D ===
plt.figure()
for nodes, values in sorted(grouped_by_nodes.items()):
    # Ordina i punti per numero di core crescente
    values.sort()
    x = [cores for cores, _ in values]
    y = [time for _, time in values]
    plt.plot(x, y, marker='o', label=f"{nodes} thread")

plt.xlabel("Numero di processi")
plt.ylabel("Speedup")
plt.title(f"Speedup vs processes per total_fishes = {selected_fishes}")
plt.legend(title="Numero di threads")
plt.grid(True)
plt.tight_layout()
plt.show()
