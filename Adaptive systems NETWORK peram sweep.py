import subprocess
import numpy as np
import matplotlib.pyplot as plt

# Define parameter values and neuron labels
P_r_values = [0.1, 0.3, 0.5, 0.7, 0.9,]
neuron_labels = [ "N2", "N3"]
sweep_results = {neuron: {} for neuron in neuron_labels}

# Directory for saving/loading files
output_dir = r"C:\Users\david\OneDrive\Desktop\PYTHONS\New folder"

# Run simulations for each P_r
for p in P_r_values:
    subprocess.run(["python", r"C:\Users\david\OneDrive\Desktop\PYTHONS\New folder\adaptive systems Network project 0.4.py", str(p)])

    # Load voltage history for each neuron
    for neuron in neuron_labels:
        filename = f"{output_dir}\Vm_history_{neuron}_{p:.2f}.npy"
        try:
            sweep_results[neuron][p] = np.load(filename)
        except FileNotFoundError:
            print(f"Warning: File not found: {filename}")
            sweep_results[neuron][p] = None  # Handle missing files gracefully

# Plot each neuron's Vm over all P_r values
for neuron in neuron_labels:
    plt.figure(figsize=(10, 6))
    for p, volt_trace in sweep_results[neuron].items():
        if volt_trace is not None:  # Skip missing data
            plt.plot(volt_trace, label=f"P_r = {p:.2f}", linewidth=1)
    plt.title(f"Membrane Potential History: {neuron}")
    plt.xlabel("Time (ms)")
    plt.ylabel("Membrane Voltage (V)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()