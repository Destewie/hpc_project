import json
import numpy as np
import matplotlib
matplotlib.use("TkAgg")  # Usa il backend TkAgg per evitare problemi con Qt/xcb
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from benchmark_functions import max_spherical_function, min_spherical_function, min_rastrigin_function

def validate_json(data):
    if not isinstance(data, list):
        raise ValueError("JSON data must be a list of epochs.")

    for epoch in data:
        if not isinstance(epoch, list):
            raise ValueError("Each epoch must be a list of fish objects.")

        for fish in epoch:
            if "x" not in fish or "weight" not in fish:
                raise ValueError("Each fish object must have 'x' and 'weight' keys.")
            if not isinstance(fish["x"], list) or len(fish["x"]) != 2:
                raise ValueError("'x' must be a list with exactly two values.")

def read_json(filepath):
    with open(filepath, 'r') as file:
        data = json.load(file)
    validate_json(data)
    return data

def create_animation(data):
    # function = min_rastrigin_function
    function = min_spherical_function

    # Determina i limiti di spawn dai dati
    all_positions = [coord for epoch in data for fish in epoch for coord in fish["x"]]
    spawn_bounds = [min(all_positions), max(all_positions)]

    # Estrarre le posizioni e i pesi per ogni epoca
    positions = [[(fish["x"][0], fish["x"][1]) for fish in epoch] for epoch in data]
    weights = [[fish["weight"] for fish in epoch] for epoch in data]

    # Trova i limiti dei pesi globali
    all_weights = [weight for epoch in weights for weight in epoch]
    min_weight, max_weight = min(all_weights), max(all_weights)

    # Configurazione del grafico
    fig, ax = plt.subplots()
    ax.set_xlim(spawn_bounds[0], spawn_bounds[1])
    ax.set_ylim(spawn_bounds[0], spawn_bounds[1])
    ax.set_title("Fish School Optimization - 2D Function")

    # Crea la griglia per la funzione
    x = np.linspace(spawn_bounds[0], spawn_bounds[1], 100)
    y = np.linspace(spawn_bounds[0], spawn_bounds[1], 100)
    X, Y = np.meshgrid(x, y)
    Z = function(X, Y)

    # Disegna la funzione
    ax.contour(X, Y, Z, 20, cmap='viridis')

    # Aggiungi i punti
    scat = ax.scatter([], [], color='red', zorder=5)

    def update(frame):
        current_positions = positions[frame]
        current_weights = weights[frame]

        # Normalizza i pesi per avere dimensioni proporzionate
        normalized_weights = [
            5 + 60 * (w - min_weight) / (max_weight - min_weight)  # Imposta una dimensione minima di 10
            for w in current_weights
        ]

        ax.clear()
        ax.set_xlim(spawn_bounds[0], spawn_bounds[1])
        ax.set_ylim(spawn_bounds[0], spawn_bounds[1])
        ax.set_title(f"Epoch {frame}")

        # Ridisegna la funzione
        ax.contour(X, Y, Z, 20, cmap='viridis')

        # Disegna i punti (i pesci)
        x_data, y_data = zip(*current_positions)
        scat = ax.scatter(x_data, y_data, s=normalized_weights, color='red', zorder=5)
        return scat,

    ani = FuncAnimation(fig, update, frames=len(positions), interval=500)
    plt.show()


# Main Execution
if __name__ == "__main__":
    try:
        # Percorso del file JSON
        filepath = "../evolution_logs/min_sphere_2d_log.json"
        # filepath = "../evolution_logs/min_rastrigin_2d_log.json"
        data = read_json(filepath)
        create_animation(data)
    except Exception as e:
        print(f"Error: {e}")
