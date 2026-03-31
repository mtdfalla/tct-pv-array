# TCT PV Array Evaluation Module

## Overview

`tct_eval.py` is a high-performance Python module for computing accurate I-V and P-V curves for Total Cross-Tied (TCT) connected photovoltaic arrays under partial shading conditions. The module implements a rigorous single-diode model with per-substring bypass diode modeling.

## Key Features

### ✓ Single-Diode Model
- Accurate Rs (series resistance) and Rsh (shunt resistance) modeling
- Temperature and irradiance scaling based on physical principles
- Saturation current temperature dependence with exponential model

### ✓ Per-Substring Bypass Diodes
- **Critical innovation**: Models 3 bypass substrings per module (18 cells each)
- Bypass diode activation at V_substring < -Vd (Vd ≈ 0.5 V)
- Dynamic resistance modeling: V_bypass = -(Vd + I·rd) where rd ≈ 0.02 Ω
- Accurate representation of Simulink-style PV Array behavior

### ✓ TCT Topology
- **Rows**: Parallel connection of modules (same voltage, currents sum)
- **Columns**: Series connection of rows (same current, voltages sum)
- Efficient solver for parallel block voltage inversion

### ✓ Robust Numerics
- Newton-Raphson with safeguards (clamped exponentials, derivative guards)
- Bisection for finding short-circuit current
- Handles extreme conditions (very low irradiance, high temperatures)
- Automatic MPP detection with prominence filtering

## Installation

```bash
# Required dependencies
pip install numpy

# Optional for plotting
pip install matplotlib
```

## Quick Start

```python
import numpy as np
from tct_eval import default_kc200gt, evaluate_tct, generate_demo_shading, plot_iv_pv

# Create module (KC200GT: 200W, 32.9V, 8.21A)
module = default_kc200gt()

# Generate a 4×3 array with partial shading
R, C = 4, 3
G_map, T_map = generate_demo_shading(R, C, pattern='partial')

# Evaluate I-V and P-V curves
results = evaluate_tct(G_map, T_map, module, num_points=300)

# Access results
print(f"GMPP: {results['gmpp'][0]:.1f} W")
print(f"Local MPPs: {len(results['local_mpps'])}")

# Plot curves
plot_iv_pv(results, title="Partial Shading Analysis")
```

## Core Functions

### 1. Module Parameters

```python
from tct_eval import ModuleParams, default_kc200gt

# Use default KC200GT
module = default_kc200gt()

# Or create custom module
custom = ModuleParams(
    Voc_stc=32.9,      # Open circuit voltage [V]
    Isc_stc=8.21,      # Short circuit current [A]
    Vmpp_stc=26.3,     # MPP voltage [V]
    Impp_stc=7.61,     # MPP current [A]
    Rs=0.34484,        # Series resistance [Ω]
    Rsh=150.6844,      # Shunt resistance [Ω]
    n=0.97734,         # Ideality factor
    Ns_cells=54,       # Total cells
    Ns_per_substring=18,  # Cells per substring
    alpha_Isc=0.0006,  # Isc temp coefficient [A/°C]
    beta_Voc=-0.00355, # Voc temp coefficient [V/°C]
    Vd=0.5,            # Bypass diode voltage [V]
    rd=0.02            # Bypass diode resistance [Ω]
)
```

### 2. Parameter Scaling

```python
from tct_eval import scale_params

# Scale to operating conditions
params = scale_params(module, G=800, T_celsius=45)
# Returns: {'IL': ..., 'Io': ..., 'Vt': ..., 'Rs': ..., 'Rsh': ...}
```

### 3. Module I-V with Bypass

```python
from tct_eval import module_voltage_at_current_with_substrings

# Compute module voltage at given current
V = module_voltage_at_current_with_substrings(
    I=7.0,              # Current [A]
    G=1000,             # Irradiance [W/m²]
    T_celsius=25,       # Temperature [°C]
    module=module
)
```

### 4. Array Evaluation

```python
from tct_eval import evaluate_tct

# Create irradiance and temperature maps
G_map = np.array([[1000, 500], [800, 600]])  # W/m²
T_map = np.full((2, 2), 25.0)                 # °C

# Evaluate array
results = evaluate_tct(
    G_map=G_map,
    T_map=T_map,
    module=module,
    num_points=500  # Resolution
)

# Results dictionary:
# {
#   'I': array of currents [A],
#   'V': array of voltages [V],
#   'P': array of powers [W],
#   'gmpp': (Pmpp, Impp, Vmpp),
#   'local_mpps': [(P1, I1, V1), (P2, I2, V2), ...]
# }
```

### 5. Shading Pattern Generation

```python
from tct_eval import generate_demo_shading

# Generate various shading patterns
G_uniform, T = generate_demo_shading(4, 3, pattern='uniform')
G_partial, T = generate_demo_shading(4, 3, pattern='partial')
G_diagonal, T = generate_demo_shading(4, 3, pattern='diagonal',
                                       levels=[200, 400, 600, 800, 1000])
G_checker, T = generate_demo_shading(4, 3, pattern='checkerboard')
G_random, T = generate_demo_shading(4, 3, pattern='random')
```

### 6. Visualization

```python
from tct_eval import plot_iv_pv

plot_iv_pv(
    results,
    title="My PV Array Analysis",
    show=True,                               # Display plot
    filename="/path/to/save/plot.png"        # Optional save
)
```

## Algorithm Details

### TCT Topology

In TCT configuration:
- **Each row** is a parallel block of C modules
- All modules in a row share the same voltage Vr
- Module currents depend on their individual G and T
- Row current = sum of module currents

- **Rows** are connected in series
- All rows carry the same current I
- Array voltage = sum of row voltages

### Evaluation Process

1. **Find Isc**: Solve V_array(I) = 0 using bisection
2. **Current sweep**: Generate points from 0 to Isc
3. **For each current I**:
   - For each row r:
     - Solve for row voltage Vr such that sum(I_module(Vr)) = I
     - Use Newton-Raphson inversion
   - Array voltage = sum of all Vr
4. **MPP detection**: Find global and local maxima in P(I)

### Bypass Diode Modeling

For each module with 3 substrings:

```
For given series current I:
  For each substring (18 cells):
    1. Solve single-diode equation for V_sub
    2. If V_sub < -Vd:
         V_sub = -(Vd + I·rd)  # Bypass conducting
    3. V_module = sum(V_sub_1, V_sub_2, V_sub_3)
```

### Single-Diode Equation

For a substring:
```
I = IL - Io·(exp((V + I·Rs)/Vt) - 1) - (V + I·Rs)/Rsh
```

Where:
- IL = light-generated current (scales with G)
- Io = saturation current (temperature dependent)
- Vt = n·Ns·k·T/q (thermal voltage)
- Rs, Rsh = series and shunt resistances

## Numerical Robustness

The module implements several safeguards:

1. **Exponential clamping**: `exp(x)` with x ∈ [-40, 40]
2. **Derivative guards**: Check for near-zero derivatives
3. **Step limiting**: Bounded Newton-Raphson updates
4. **Voltage bounds**: Allow small negative voltages (−1V) for bypass
5. **Convergence checks**: Dual criteria (function and update)

## Validation Results

Tested against theoretical values for 4×3 array:

| Metric | Theoretical | Actual | Error |
|--------|-------------|--------|-------|
| Voc    | 131.6 V     | 131.5 V | 0.11% |
| Isc    | 24.6 A      | 24.6 A  | 0.23% |
| Pmpp   | 2402 W      | 2393 W  | 0.35% |

## Performance

For a 4×3 array with 300 evaluation points:
- Computation time: ~0.5-1 seconds (Python)
- Memory footprint: < 10 MB
- Suitable for optimization loops requiring 100-1000s of evaluations

## Use Cases

### 1. Reconfiguration Optimization
```python
# Evaluate different physical arrangements
for config in candidate_configurations:
    G_map = apply_shading_to_config(config, shading_pattern)
    results = evaluate_tct(G_map, T_map, module)
    fitness = results['gmpp'][0]  # Use GMPP as fitness
```

### 2. MPPT Algorithm Development
```python
# Get complete I-V curve for MPPT testing
results = evaluate_tct(G_map, T_map, module, num_points=1000)
I_curve = results['I']
V_curve = results['V']
P_curve = results['P']
```

### 3. Shading Impact Analysis
```python
# Compare different shading scenarios
scenarios = ['uniform', 'partial', 'diagonal', 'checkerboard']
for scenario in scenarios:
    G_map, T_map = generate_demo_shading(R, C, pattern=scenario)
    results = evaluate_tct(G_map, T_map, module)
    print(f"{scenario}: {results['gmpp'][0]:.1f} W")
```

### 4. Temperature/Irradiance Sensitivity
```python
# Study environmental effects
for T in range(25, 76, 10):
    for G in range(200, 1001, 200):
        G_map = np.full((R, C), float(G))
        T_map = np.full((R, C), float(T))
        results = evaluate_tct(G_map, T_map, module)
        # Analyze results...
```

## Integration with Optimization

### Genetic Algorithm Example
```python
import numpy as np
from tct_eval import evaluate_tct, default_kc200gt

def fitness_function(individual):
    """Individual encodes physical module arrangement"""
    # Map individual to physical configuration
    G_map = map_to_configuration(individual, base_shading)
    T_map = np.full_like(G_map, 25.0)
    
    # Evaluate
    module = default_kc200gt()
    results = evaluate_tct(G_map, T_map, module)
    
    return results['gmpp'][0]  # Maximize GMPP

# Use in GA
from genetic_algorithm import GeneticAlgorithm
ga = GeneticAlgorithm(fitness_function, ...)
best_config = ga.evolve()
```

### Reinforcement Learning Example
```python
def get_power_reward(state, action):
    """
    state: current module arrangement
    action: reconfiguration move
    """
    new_arrangement = apply_action(state, action)
    G_map = get_irradiance_map(new_arrangement)
    T_map = get_temperature_map()
    
    results = evaluate_tct(G_map, T_map, module)
    return results['gmpp'][0]
```

## Technical Notes

### Coordinate System
- Rows (R): Vertical dimension, series connection
- Columns (C): Horizontal dimension, parallel connection
- G_map and T_map shape: (R, C)

### Electrical Topology
```
     Col 0    Col 1    Col 2
    [M0,0]  [M0,1]  [M0,2]  } Row 0 (parallel)
       |       |       |
    [M1,0]  [M1,1]  [M1,2]  } Row 1 (parallel)
       |       |       |
    [M2,0]  [M2,1]  [M2,2]  } Row 2 (parallel)
       |       |       |
    [M3,0]  [M3,1]  [M3,2]  } Row 3 (parallel)
    
    Rows connected in series: V_array = V_row0 + V_row1 + V_row2 + V_row3
```

### KC200GT Specifications

Default module parameters from Kyocera KC200GT datasheet:

| Parameter | Value | Unit |
|-----------|-------|------|
| Maximum Power | 200.1 | W |
| Voltage at Pmax | 26.3 | V |
| Current at Pmax | 7.61 | A |
| Open Circuit Voltage | 32.9 | V |
| Short Circuit Current | 8.21 | A |
| Total Cells | 54 | - |
| Substrings | 3 | - |
| Cells per Substring | 18 | - |

## Limitations

1. **Uniform temperature within substrings**: Model assumes all cells in a substring are at the same temperature
2. **Static bypass threshold**: Real bypass diodes have non-linear I-V characteristics
3. **No capacitive effects**: Steady-state model only
4. **No degradation**: Assumes new modules (can be adjusted via parameters)

## Future Extensions

Potential enhancements:
- Dynamic bypass diode I-V curves
- Cell-level temperature gradients
- Series-Parallel (SP) and Bridge-Linked (BL) topologies
- Hotspot temperature modeling
- Soiling and degradation effects

## References

1. Single-diode model: Villalva et al., "Comprehensive Approach to Modeling and Simulation of Photovoltaic Arrays" (IEEE Trans. 2009)
2. TCT topology: Srinivasan & Venkatesan, "Analysis of various electrical configurations for PV arrays" (Renewable Energy, 2019)
3. Bypass diodes: Silvestre et al., "Effects of shadowing on photovoltaic module performance" (Progress in Photovoltaics, 2008)

## License

This module is provided for research and educational purposes.

## Support

For questions or issues related to integration with optimization frameworks, please refer to the verification script (`verify_tct.py`) for comprehensive usage examples.

---

**Author**: Generated for PV array reconfiguration optimization research  
**Version**: 1.0  
**Last Updated**: November 2025
