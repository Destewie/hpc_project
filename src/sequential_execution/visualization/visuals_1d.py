import json
import numpy as np
import matplotlib
matplotlib.use("TkAgg")  # Usa il backend TkAgg per evitare problemi con Qt/xcb
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from benchmark_functions import min_spherical_function, max_spherical_function, min_rastringin_function

def validate_json(data):
    if not isinstance(data, list):
        raise ValueError("JSON data must be a list of epochs.")

    for epoch in data:
        if not isinstance(epoch, list):
            raise ValueError("Each epoch must be a list of fish objects.")

        for fish in epoch:
            if "x" not in fish or "weight" not in fish:
                raise ValueError("Each fish object must have 'x' and 'weight' keys.")
            if not isinstance(fish["x"], list) or len(fish["x"]) != 1:
                raise ValueError("'x' must be a list with a single value.")

def read_json(filepath):
    with open(filepath, 'r') as file:
        data = json.load(file)
    validate_json(data)
    return data

def create_animation(data):
    # function = max_spherical_function
    function = min_rastringin_function

    # Determina i limiti di spawn dai dati
    all_positions = [fish["x"][0] for epoch in data for fish in epoch]
    spawn_bounds = [min(all_positions), max(all_positions)]

    # Estrarre le posizioni e i pesi per ogni epoca
    positions = [[fish["x"][0] for fish in epoch] for epoch in data]
    weights = [[fish["weight"] for fish in epoch] for epoch in data]

    # Configurazione del grafico
    fig, ax = plt.subplots()
    ax.set_xlim(spawn_bounds[0], spawn_bounds[1])
    ax.set_ylim(0, function(spawn_bounds[1]))  # Imposta i limiti in y in base alla funzione
    ax.set_title("Fish School Optimization - Unidimensional Function")

    # Crea il vettore di valori x per il grafico della funzione
    x = np.linspace(spawn_bounds[0], spawn_bounds[1], 100)
    y = function(x)

    # Disegna la funzione
    ax.plot(x, y, label="Fitness Function", color='blue')

    # Aggiungi i punti
    scat = ax.scatter([], [], color='red', zorder=5)

    def update(frame):
        current_positions = positions[frame]
        current_weights = weights[frame]

        # Normalizza i pesi per avere dimensioni proporzionate
        normalized_weights = np.array(current_weights) * 20

        ax.clear()
        ax.set_xlim(spawn_bounds[0], spawn_bounds[1])
        ax.set_ylim(0, function(spawn_bounds[1]))  # Imposta i limiti in y in base alla funzione
        ax.set_title(f"Epoch {frame}")

        # Ridisegna la funzione
        ax.plot(x, y, label="Fitness Function", color='blue')

        # Disegna i punti (i pesci)
        y_data = function(np.array(current_positions))  # Calcola i valori di y per i punti
        scat = ax.scatter(current_positions, y_data, s=normalized_weights, color='red', zorder=5)
        return scat,

    ani = FuncAnimation(fig, update, frames=len(positions), interval=500)
    plt.show()

# Main Execution
if __name__ == "__main__":
    try:
        # Percorso del file JSON
        filepath = "../evolution_logs/min_rastringin_1d_log.json"
        data = read_json(filepath)
        create_animation(data)
    except Exception as e:
        print(f"Error: {e}")
