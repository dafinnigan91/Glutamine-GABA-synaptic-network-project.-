# Glutamine-GABA-synaptic-network-project.-

# Synaptic Dynamics Simulation

A comprehensive computational model of glutamatergic and GABAergic synaptic transmission with dynamic ion concentrations, receptor kinetics, and synaptic plasticity mechanisms.

## Overview

This project implements a detailed biophysical simulation of synaptic function, modeling the complex interactions between neurotransmitter release, receptor binding, ion flux, and membrane potential dynamics. The simulation captures both excitatory (glutamate) and inhibitory (GABA) synaptic transmission with realistic plasticity mechanisms.

## Features

### Biophysical Accuracy
- **Goldman-Hodgkin-Katz equation** for dynamic membrane potential calculation
- ** ion concentrations** (Na⁺, K⁺, Cl⁻, Ca²⁺) with proper equilibrium values
- **Vesicle release dynamics** with probabilistic neurotransmitter release
- **Receptor kinetics** for AMPA, NMDA, GABA-A, and GABA-B receptors

### Synaptic Transmission
- **Glutamatergic synapses**: AMPAR and NMDAR with Mg²⁺ block mechanism
- **GABAergic synapses**: GABA-A and GABA-B receptor dynamics
- **Calcium-dependent plasticity** (LTP/LTD mechanisms)
- **Activity-dependent receptor trafficking**

### Physiological Mechanisms
- **Ion pumps** for membrane potential stabilization
- **Passive ion decay** and leak currents
- **Neurotransmitter clearance** by transporters
- **Synaptic depression** and facilitation

## Scientific Background

### Glutamatergic Synapses
The simulation models excitatory synaptic transmission through:
- **AMPA receptors**: Fast excitatory currents (Na⁺/K⁺)
- **NMDA receptors**: Slow, Ca²⁺-permeable currents with voltage dependence
- **Calcium-dependent plasticity**: Receptor insertion/removal based on Ca²⁺ thresholds

### GABAergic Synapses  
Inhibitory transmission is modeled via:
- **GABA-A receptors**: Fast inhibitory currents (Cl⁻)
- **GABA-B receptors**: Slow, metabotropic K⁺ currents
- **Activity-dependent inhibitory plasticity**

## Installation

### Requirements
```bash
pip install numpy matplotlib pandas
```

### Dependencies
- `numpy` - Numerical computations
- `matplotlib` - Visualization and plotting
- `pandas` - Data handling and export
- `random` - Stochastic processes
- `sys` - Parameter input handling

## Usage

### Basic Simulation
```bash
python synapse_simulation.py [vesicle_release_probability]
```

### Parameter Sweep
```bash
# Run with different release probabilities
python synapse_simulation.py 0.8
python synapse_simulation.py 0.6
python synapse_simulation.py 0.4
```

### Example Output
The simulation generates:
- **Time series plots** of membrane potential, receptor dynamics, and ion concentrations
- **Excel export** of all tracked variables
- **NumPy arrays** for further analysis

## Key Parameters

### Synaptic Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| `P_r` | 0.8 | Vesicle release probability |
| `GLu_mol_Vescl` | 3000 | Glutamate molecules per vesicle |
| `AMPAR_initial` | 500 | Initial AMPA receptor count |
| `NMDAR` | 50 | NMDA receptor count |
| `Ca_threshold` | 40 | Calcium threshold for plasticity |

### Ion Concentrations (mM)
| Ion | Intracellular | Extracellular |
|-----|---------------|---------------|
| Na⁺ | 15.0 | 145.0 |
| K⁺ | 127.0 | 5.0 |
| Cl⁻ | 4.0 | 110.0 |
| Ca²⁺ | 0.0001 | 0.1 |

## Model Architecture

### Class Structure
```python
class Gultamatergic_Synaps:
    - AMPAR/NMDAR binding and kinetics
    - Calcium-dependent plasticity
    - Ion flux calculations
    - Receptor trafficking

class GABAergic_Synaps:
    - GABA-A/GABA-B receptor dynamics
    - Inhibitory plasticity mechanisms
    - Chloride and potassium currents
```

### Key Functions
- `Membrain_potential()` - Goldman-Hodgkin-Katz calculation
- `Ion_Pump()` - Na⁺/K⁺ ATPase simulation
- `Ion_decay()` - Passive ion equilibration
- `GLu_leak()` - Neurotransmitter clearance

## Output Data

### Time Series Data
- Membrane potential (mV)
- Ion concentrations (mM)
- Receptor numbers
- Spike trains
- Synaptic currents

### Visualization
- Multi-panel time series plots
- Ion concentration dynamics
- Receptor trafficking over time
- Membrane potential evolution

## Applications

### Research Applications
- **Synaptic plasticity studies** - LTP/LTD mechanisms
- **Pharmacological modeling** - Drug effects on receptors
- **Disease modeling** - Synaptic dysfunction in neurological disorders
- **Educational tool** - Understanding synaptic physiology

### Computational Neuroscience
- **Network integration** - Use as building block for neural networks
- **Parameter exploration** - Systematic parameter space investigation
- **Model validation** - Compare with experimental data

## Future Developments

### Planned Features
- *Metabotropic receptor signaling** (mGluR, mAChR)
- *Calcium buffering systems** (calbindin, calmodulin)
- *Presynaptic plasticity** (paired-pulse facilitation/depression)
- *Multi-compartment integration** (dendritic tree modeling)
- *Stochastic channel dynamics** (single-channel noise)

### Model Extensions
- *Temperature dependence** of kinetic rates
- *pH sensitivity** of receptors and channels
- *Metabolic constraints** on ATP-dependent processes
- *Glial cell interactions** (astrocyte glutamate uptake)

## Contributing

Contributions are welcome! Areas of interest:
- **Experimental validation** with patch-clamp data
- **Parameter optimization** using genetic algorithms
- **Performance optimization** for large-scale simulations
- **Documentation improvements** and code cleanup

## References

### Key Publications
1. **Goldman-Hodgkin-Katz equation**: Hodgkin & Katz (1949) *J Physiol*
2. **NMDA receptor kinetics**: Jahr & Stevens (1990) *Nature*
3. **Synaptic plasticity**: Malenka & Bear (2004) *Neuron*
4. **GABAergic transmission**: Farrant & Nusser (2005) *Nat Rev Neurosci*

### Model Validation
- Ion concentration ranges from Kandel et al. *Principles of Neural Science*
- Receptor kinetics from experimental patch-clamp studies
- Plasticity mechanisms from hippocampal slice preparations

## License

This project is open source. Please cite if used in academic work:

```
Finnigan, D.A. (2025). Synaptic Dynamics Simulation: A Biophysical Model of 
Glutamatergic and GABAergic Transmission. GitHub Repository.
```

## Contact

**David A. Finnigan**  
MSc AI & Adaptive Systems
University of Sussex

- GitHub: [@dafinnigan91](https://github.com/dafinnigan91)
- LinkedIn: [david-finnigan-ai](https://linkedin.com/in/david-finnigan-ai)

---


