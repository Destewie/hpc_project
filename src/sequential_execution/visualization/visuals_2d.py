import json
import numpy as np
import matplotlib
matplotlib.use("TkAgg")  # Usa il backend TkAgg per evitare problemi con Qt/xcb
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Funzione per validare la struttura JSON
def validate_json(data):
    required_keys = ["fitness_function", "epochs"]
    for key in required_keys:
        if key not in data:
            raise ValueError(f"Missing key: '{key}' in JSON.")

    fitness_function_keys = ["type", "formula", "n_dimensions", "spawn_bounds"]
    for key in fitness_function_keys:
        if key not in data["fitness_function"]:
            raise ValueError(f"Missing key: '{key}' in 'fitness_function' section.")

    if not isinstance(data["fitness_function"]["n_dimensions"], int) or data["fitness_function"]["n_dimensions"] != 2:
        raise ValueError("'n_dimensions' must be 2 for this visualization.")

    if len(data["fitness_function"]["spawn_bounds"]) != 2:
        raise ValueError("'spawn_bounds' must contain exactly 2 values.")

# Funzione per leggere e validare il file JSON
def read_json(filepath):
    with open(filepath, 'r') as file:
        data = json.load(file)
    validate_json(data)
    return data

# Funzione per calcolare la funzione sferica
def spherical_function(x, y):
    return x**2 + y**2

# Funzione per creare l'animazione
def create_animation(data):
    spawn_bounds = data["fitness_function"]["spawn_bounds"]
    positions = [[(fish["x"], fish["y"]) for fish in epoch] for epoch in data["epochs"]]

    # Configurazione del grafico
    fig, ax = plt.subplots()
    ax.set_xlim(spawn_bounds[0], spawn_bounds[1])
    ax.set_ylim(spawn_bounds[0], spawn_bounds[1])
    ax.set_title("Fish School Optimization - Spherical Function")

    # Genera la griglia per la funzione
    x = np.linspace(spawn_bounds[0], spawn_bounds[1], 100)
    y = np.linspace(spawn_bounds[0], spawn_bounds[1], 100)
    X, Y = np.meshgrid(x, y)
    Z = spherical_function(X, Y)

    # Aggiungi i contorni della funzione (senza mappa di colore)
    ax.contour(X, Y, Z, levels=20, cmap="viridis", alpha=0.7)

    # Aggiungi i punti
    scat = ax.scatter([], [], color='red', zorder=5)

    def update(frame):
        current_positions = positions[frame]
        ax.clear()
        ax.set_xlim(spawn_bounds[0], spawn_bounds[1])
        ax.set_ylim(spawn_bounds[0], spawn_bounds[1])
        ax.set_title(f"Epoch {frame}")

        # Ridisegna i contorni della funzione
        ax.contour(X, Y, Z, levels=20, cmap="viridis", alpha=0.7)

        # Disegna i punti
        x_data, y_data = zip(*current_positions)
        scat = ax.scatter(x_data, y_data, color='red', zorder=5)
        return scat,

    ani = FuncAnimation(fig, update, frames=len(positions), interval=500)
    plt.show()

# Main Execution
if __name__ == "__main__":
    try:
        # Percorso del file JSON
        filepath = "../evolution_logs/spherical_2d_log.json"
        data = read_json(filepath)
        create_animation(data)
    except Exception as e:
        print(f"Error: {e}")
