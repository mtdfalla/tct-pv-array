"""
Test script for 9×9 array with custom irradiance pattern
"""

import numpy as np
from tct_eval import default_kc200gt, evaluate_tct, plot_iv_pv

print("=" * 70)
print("9×9 Array Evaluation - Custom Irradiance Pattern")
print("=" * 70)

# Your custom 9×9 irradiance matrix
G = np.array([
    [100,  600,  400,  400,  400,  400,  100,  400,  600],
    [600,  600, 1000,  800, 1000,  800, 1000,  600,  600],
    [1000, 400, 1000,  400,  600, 1000,  600,  600,  600],
    [200,  400,  100,  200,  200,  600,  100,  200,  200],
    [1000, 800, 1000,  600,  600,  800,  800,  400, 1000],
    [600,  400,  600,  800, 1000,  800,  800,  400,  800],
    [600,  600,  100,  100,  600,  100,  600,  200,  600],
    [1000, 400, 1000, 1000, 1000,  100,  800,  800, 1000],
    [1000, 100,  800,  400,  800,  200,  600,  800,  800]
])



# G = np.array([
#     [100,  600,  400,  400,  400,  400,  100,  400,  600],
#     [200,  600, 1000,  800, 1000,  800, 1000,  600,  600],
#     [900, 400, 1000,  100,  600, 1000,  100,  600,  600],
#     [200,  400,  100,  200,  200,  600,  100,  200,  200],
#     [100, 800, 1000,  600,  600,  800,  800,  400, 1000],
#     [600,  400,  600,  800, 1000,  800,  800,  400,  800],
#     [600,  600,  100,  100,  600,  100,  600,  200,  600],
#     [300, 400, 1000, 1000, 1000,  100,  800,  800, 1000],
#     [1000, 100,  800,  0,  800,  200,  600,  800,  800]
# ])



# G = np.array([
#     [1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000],
#     [1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000],
#     [1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000],
#     [1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000],
#     [1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000],
#     [1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000],
#     [1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000],
#     [1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000],
#     [1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000]
# ])

# Temperature (uniform at 25°C)
T = np.full((9, 9), 25.0)

# Display the irradiance pattern
print("\nIrradiance Map [W/m²]:")
print(G)
print(f"\nArray size: {G.shape[0]}×{G.shape[1]} = {G.size} modules")
print(f"Average irradiance: {G.mean():.1f} W/m²")
print(f"Min: {G.min()} W/m², Max: {G.max()} W/m²")

# Load module parameters
module = default_kc200gt()
print(f"\nModule: KC200GT (200W per module at STC)")
print(f"Total rated power (STC): {G.size * 200:.0f} W")

# Evaluate the array
print("\nEvaluating I-V and P-V curves...")
results = evaluate_tct(G, T, module, num_points=500)

# Display results
print("\n" + "=" * 70)
print("RESULTS")
print("=" * 70)

gmpp = results['gmpp']
print(f"\nGlobal Maximum Power Point (GMPP):")
print(f"  Power:   {gmpp[0]:.2f} W")
print(f"  Current: {gmpp[1]:.2f} A")
print(f"  Voltage: {gmpp[2]:.2f} V")

print(f"\nArray Characteristics:")
print(f"  Open Circuit Voltage (Voc): {results['V'][0]:.2f} V")
print(f"  Short Circuit Current (Isc): {results['I'][-1]:.2f} A")

print(f"\nLocal Maximum Power Points:")
print(f"  Total local MPPs detected: {len(results['local_mpps'])}")
print(f"\nTop 5 Power Peaks:")
for i, (P, I, V) in enumerate(results['local_mpps'][:5]):
    print(f"  {i+1}. P = {P:.2f} W, I = {I:.2f} A, V = {V:.2f} V")

# Calculate efficiency under partial shading
ideal_power = G.size * 200  # All modules at STC
actual_power = gmpp[0]
efficiency = (actual_power / ideal_power) * 100
print(f"\nArray Performance:")
print(f"  Ideal power (all modules at STC): {ideal_power:.0f} W")
print(f"  Actual GMPP: {actual_power:.2f} W")
print(f"  Efficiency under shading: {efficiency:.1f}%")

# Generate and save plots
print("\n" + "=" * 70)
print("Generating plots...")
print("=" * 70)

plot_iv_pv(
    results, 
    title="9×9 Array - Custom Irradiance Pattern",
    filename="9x9_custom_pattern.png",
    show=False
)

print("\nPlot saved as: 9x9_custom_pattern.png")
print("\n" + "=" * 70)
print("Evaluation Complete!")
print("=" * 70)