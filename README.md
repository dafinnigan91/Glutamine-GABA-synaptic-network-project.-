# Neural Network Synaptic Modeling: Chaos and Adaptation in Silico

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/release/python-3110/)
[![NumPy](https://img.shields.io/badge/numpy-%23013243.svg?style=flat&logo=numpy&logoColor=white)](https://numpy.org/)
[![NetworkX](https://img.shields.io/badge/networkx-2.0+-orange.svg)](https://networkx.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 🧠 Overview

A computational model of **biologically-inspired neural networks** featuring both excitatory (glutamatergic) and inhibitory (GABAergic) synapses. This project explores **adaptive behavior** and **chaotic dynamics** in artificial neural circuits through parameter analysis and Lyapunov exponent characterization.

### Key Features
- 🔬 **Biologically accurate** synaptic mechanisms (AMPA/NMDA/GABA-A/GABA-B receptors)
- ⚡ **Dynamic plasticity** with calcium-dependent LTP/LTD
- 🌐 **Network-level coordination** using NetworkX graph structures
- 📊 **Chaos analysis** via Lyapunov exponent calculation
- 🎛️ **Parameter sweeps** for vesicle release probability analysis
- ⚖️ **Membrane potential** computation using Goldman-Hodgkin-Katz equation

## 🎯 Research Significance

This work bridges **computational neuroscience**, **network science**, and **adaptive systems theory** to understand how biological coordination mechanisms can inform artificial neural architectures.

### Applications
- **Neuromorphic AI** system design
- **Bio-inspired robotics** control systems
- **Adaptive learning** algorithm development
- **Neural network** stability analysis

## 🏗️ System Architecture

### Network Structure
```
N1 (Input) ──[Glu]──> N2 (GABAergic) inhibition
     │                    │
     └───[Glu]──> N3 ──[Glu]──> N4 (Output)
                           
                 
```

### Synaptic Components
- **Glutamatergic Synapses**: AMPA/NMDA receptor dynamics with calcium-dependent plasticity
- **GABAergic Synapses**: GABA-A/GABA-B mediated inhibition with homeostatic regulation
- **Ion Dynamics**: Na⁺, K⁺, Cl⁻, Ca²⁺ concentration modeling with realistic pumps and leaks

## 📈 Key Results

### Adaptive Behavior
- ✅ **Synaptic plasticity** responds to calcium thresholds
- ✅ **Homeostatic regulation** prevents runaway dynamics  
- ✅ **Noise robustness** integrates perturbations appropriately
- ✅ **Realistic inhibition** demonstrates IPSP suppression

### Chaos Analysis
- **Lyapunov exponents** reveal sensitivity to vesicle release probability
- **Parameter sweeps** show transitions between chaotic and stable regimes
- **Lower release probabilities** → Higher sensitivity (λ ≈ 8.208)
- **Higher release probabilities** → Increased stability

## 🚀 Quick Start

### Prerequisites
```bash
pip install numpy matplotlib networkx pandas
```

### Basic Usage
```python
from neural_network_model import *

# Initialize network with 4 neurons
network = create_neural_network()

# Run simulation with default parameters
results = run_simulation(timesteps=9999, release_prob=0.8)

# Analyze results
plot_membrane_potentials(results)
plot_receptor_dynamics(results)
calculate_lyapunov_exponent(results)
```

### Parameter Sweep Analysis
```python
# Test range of release probabilities
prob_range = np.linspace(0.01, 1.0, 20)
lyapunov_results = []

for prob in prob_range:
    lyap = compute_lyapunov_exponent(prob)
    lyapunov_results.append(lyap)

plot_chaos_analysis(prob_range, lyapunov_results)
```

## 📁 Project Structure

```
neural-synaptic-modeling/
├── src/
│   ├── synaptic_models.py          # Core synapse classes
│   ├── network_simulation.py      # NetworkX integration
│   ├── analysis_tools.py           # Lyapunov & parameter analysis
│   └── visualization.py           # Plotting functions
├── data/
│   ├── baseline_parameters.json   # Default ion concentrations
│   └── results/                   # Simulation outputs
├── docs/
│   ├── mathematical_framework.md  # GHK equation derivations
│   ├── biological_background.md   # Synaptic mechanisms
│   └── chaos_theory.md           # Lyapunov analysis theory
├── tests/
│   ├── test_synapses.py           # Unit tests
│   └── test_network.py            # Integration tests
└── examples/
    ├── basic_simulation.py        # Getting started
    ├── parameter_sweep.py         # Chaos analysis
    └── custom_network.py          # Advanced usage
```

## 🧮 Mathematical Framework

### Goldman-Hodgkin-Katz Equation
```
Vm = (RT/F) × ln[(PK[K+]out + PNa[Na+]out + PCl[Cl-]in + PCa[Ca2+]out) / 
                  (PK[K+]in + PNa[Na+]in + PCl[Cl-]out + PCa[Ca2+]in)]
```

### Lyapunov Exponent Calculation
```python
def compute_lyapunov_exponent(binding_prob, delta=1e-5, timesteps=9999):
    """
    Quantifies system sensitivity to initial conditions
    
    Returns:
        float: Lyapunov exponent (λ > 0 indicates chaos)
    """
    trajectory1 = simulate_system(binding_prob)
    trajectory2 = simulate_system(binding_prob + delta)
    
    divergence = np.abs(trajectory1 - trajectory2)
    divergence[divergence == 0] = 1e-8
    
    return np.mean(np.log(divergence / delta))
```

## 🔬 Research Methods

### Synaptic Plasticity Rules
- **LTP Induction**: Ca²⁺ > threshold → AMPAR insertion (200-300 units)
- **LTD Induction**: Prolonged low Ca²⁺ → AMPAR removal (50-100 units) 
- **Homeostatic Bounds**: Receptor counts clamped to realistic ranges
- **Vesicle Dynamics**: Stochastic release with fatigue and recovery

### Network Validation
- **Signal Propagation**: Realistic excitatory transmission N1→N2, N1→N3
- **Inhibitory Control**: GABAergic suppression N2→N4 prevents runaway excitation
- **Membrane Stability**: Voltage dynamics remain within physiological bounds
- **Noise Integration**: Random perturbations absorbed without destabilization

## 📊 Analysis Tools

### Visualization Functions
```python
# Membrane potential traces
plot_membrane_dynamics(neuron_data, neuron_id='N2')

# Receptor evolution over time  
plot_receptor_plasticity(synapse_data, receptor_type='AMPAR')

# Parameter sensitivity analysis
plot_parameter_sweep(prob_range, membrane_responses)

# Chaos characterization
plot_lyapunov_analysis(probabilities, exponents)
```

### Data Export
- **NumPy arrays** for numerical analysis
- **JSON format** for parameter configurations  
- **CSV export** for external statistical analysis
- **Figure generation** for publication-ready plots

## 🎓 Educational Value

### Learning Objectives
- **Computational Neuroscience**: Synaptic mechanisms and plasticity
- **Network Science**: Graph-based neural connectivity
- **Chaos Theory**: Sensitivity analysis and stability transitions  
- **Biological Modeling**: Ion dynamics and membrane potentials
- **Python Programming**: Scientific computing and visualization

### Extensions for Students
1. **Add new receptor types** (mGluR, GABAB metabotropic)
2. **Implement action potentials** (Hodgkin-Huxley dynamics)
3. **Expand network size** (scale to 100+ neurons)
4. **Add spatial dynamics** (dendritic compartments)
5. **Include neuromodulation** (dopamine, serotonin effects)

## 🤝 Contributing

We welcome contributions! Areas of particular interest:

- **Enhanced biological realism** (SNARE proteins, G-protein cascades)
- **Performance optimization** (vectorized computations, GPU acceleration)
- **Additional analysis tools** (mutual information, transfer entropy)
- **Visualization improvements** (3D network plots, interactive dashboards)

### Development Setup
```bash
git clone https://github.com/dafinnigan91/neural-synaptic-modeling.git
cd neural-synaptic-modeling
pip install -r requirements.txt
pip install -e .  # Install in development mode

# Run tests
python -m pytest tests/
```

## 📚 References & Background

### Key Publications
- **Synaptic Plasticity**: Hebbian learning and calcium-dependent mechanisms
- **Chaos in Neural Networks**: Lyapunov analysis of biological circuits  
- **Network Neuroscience**: Graph theory applications to brain connectivity
- **Computational Models**: Biologically-inspired artificial neural systems

### Research Context
This project contributes to understanding how **biological coordination mechanisms** can inform **artificial intelligence architectures**. The chaos analysis reveals how **release probability** acts as a critical control parameter governing network stability—insights applicable to multi-agent AI systems and neuromorphic computing.

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 👨‍🔬 Author

**David Finnigan**  
MSc AI & Adaptive Systems (Distinction), University of Sussex

- 🔗 Research Focus: Computational Neuroscience, Network Science, Multi-Agent Coordination
- 💻 GitHub: [@dafinnigan91](https://github.com/dafinnigan91)
- 📧 Contact: [Your Email]

---

## 🌟 Acknowledgments

- **University of Sussex** Adaptive Systems Program
- **NetworkX Community** for graph-based modeling tools
- **Computational Neuroscience** research community
- **Open Source Scientific Python** ecosystem

---

*This project bridges biological neural mechanisms with artificial intelligence, exploring how chaos and adaptation in synaptic networks can inform next-generation AI architectures.*
