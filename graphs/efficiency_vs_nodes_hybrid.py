import json
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from collections import defaultdict

# === Leggi il file JSON ===
with open("../implementations/hybrid/results_modified.json", "r") as f:
    data = json.load(f)

# === Trova i valori unici di total_fishes ===
fishes_options = sorted({entry["total_fishes"] for entry in data.values()})
print("Valori disponibili per 'total_fishes':")
for i, val in enumerate(fishes_options):
    print(f"[{i}] {val}")

# === Scelta total_fishes ===
try:
    index = int(input("Seleziona l'indice del valore di 'total_fishes' da visualizzare: "))
    selected_fishes = fishes_options[index]
except (IndexError, ValueError):
    print("Indice non valido.")
    exit()

# === Trova i valori unici di places per quel numero di pesci ===
places_options = sorted({entry["places"] for entry in data.values() if entry["total_fishes"] == selected_fishes})
print("\nValori disponibili per 'places':")
for i, val in enumerate(places_options):
    print(f"[{i}] {val}")

# === Scelta places ===
try:
    places_index = int(input("Seleziona l'indice del placement ('places') da visualizzare: "))
    selected_places = places_options[places_index]
except (IndexError, ValueError):
    print("Indice non valido.")
    exit()

# === Filtra i dati ===
filtered = [entry for entry in data.values()
            if entry["total_fishes"] == selected_fishes and entry["places"] == selected_places]

# === Raggruppa per numero di thread (cores per nodo) ===
grouped_by_threads = defaultdict(list)
for entry in filtered:
    threads = entry["cores"]
    grouped_by_threads[threads].append((entry["nodes"], entry["efficiency"]))

# === Crea grafico ===
plt.figure()
for threads, values in sorted(grouped_by_threads.items()):
    if threads in [1,2,4,8,16,32,64]:
        values.sort()
        x = [nodes for nodes, _ in values]
        y = [eff for _, eff in values]
        plt.plot(x, y, marker='o', label=f"{threads} thread")

plt.xticks([1, 2, 4, 8, 16, 32, 64])
plt.xlim(0, 67)
plt.ylim(0, 1.3)  # Imposta un limite superiore per l'efficienza
plt.xlabel("Numero di processi (nodes)")
plt.ylabel("Efficiency")
plt.title(f"Efficiency vs processes - fishes={selected_fishes}, places={selected_places}")
plt.legend(title="Numero di threads")
plt.grid(True)
plt.tight_layout()
# plt.show()
plt.savefig(f"images/efficiency_vs_processes_hybrid_strong_{selected_fishes}_{selected_places}.png")
