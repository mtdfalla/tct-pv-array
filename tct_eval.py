"""
TCT (Total Cross-Tied) PV Array Evaluation Module

This module provides fast, accurate I-V and P-V curve computation for TCT-connected
PV arrays under partial shading conditions. It uses a single-diode model with 
per-substring bypass diode modeling (3 substrings of 18 cells each per module).

Key features:
- Single-diode model with Rs, Rsh, and temperature/irradiance scaling
- Per-substring bypass diode activation (Vd ≈ 0.5 V, rd ≈ 0.02 Ω)
- Robust numerical methods (Newton-Raphson with safeguards)
- Efficient vectorized operations where possible
- TCT topology: rows are parallel blocks, rows connected in series
"""

import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
import warnings

# Physical constants
K_BOLTZMANN = 1.380649e-23  # J/K
Q_ELECTRON = 1.602176634e-19  # C

@dataclass
class ModuleParams:
    """Parameters for a PV module using single-diode model"""
    # Electrical parameters at STC (1000 W/m², 25°C)
    Voc_stc: float      # Open circuit voltage [V]
    Isc_stc: float      # Short circuit current [A]
    Vmpp_stc: float     # Voltage at MPP [V]
    Impp_stc: float     # Current at MPP [A]
    
    # Single-diode model parameters
    Rs: float           # Series resistance [Ω]
    Rsh: float          # Shunt resistance [Ω]
    n: float            # Ideality factor (diode quality factor)
    
    # Physical parameters
    Ns_cells: int       # Number of cells in series
    Ns_per_substring: int = 18  # Cells per substring (for bypass modeling)
    
    # Temperature coefficients
    alpha_Isc: float = 0.0006    # Temperature coefficient of Isc [A/°C]
    beta_Voc: float = -0.00355   # Temperature coefficient of Voc [V/°C]
    
    # Bypass diode parameters
    Vd: float = 0.5      # Bypass diode forward voltage [V]
    rd: float = 0.02     # Bypass diode dynamic resistance [Ω]
    
    # Reference conditions
    G_stc: float = 1000.0   # Standard irradiance [W/m²]
    T_stc: float = 25.0     # Standard temperature [°C]


def default_kc200gt() -> ModuleParams:
    """
    Return default parameters for Kyocera KC200GT module.
    
    Specifications at STC (1000 W/m², 25°C, AM 1.5):
    - Voc = 32.9 V
    - Isc = 8.21 A
    - Vmpp = 26.3 V
    - Impp = 7.61 A
    - Pmpp = 200.1 W
    - 54 cells (3 substrings × 18 cells)
    """
    return ModuleParams(
        Voc_stc=32.9,
        Isc_stc=8.21,
        Vmpp_stc=26.3,
        Impp_stc=7.61,
        Rs=0.34484,
        Rsh=150.6844,
        n=0.97734,
        Ns_cells=54,
        Ns_per_substring=18,
        alpha_Isc=0.0006,
        beta_Voc=-0.00355,
        Vd=0.5,
        rd=0.02
    )


def scale_params(module: ModuleParams, G: float, T_celsius: float) -> Dict[str, float]:
    """
    Scale module parameters from STC to operating conditions (G, T).
    
    Args:
        module: Module parameters at STC
        G: Irradiance [W/m²]
        T_celsius: Cell temperature [°C]
    
    Returns:
        Dictionary with scaled parameters:
        - IL: Light-generated current [A]
        - Io: Diode saturation current [A]
        - Vt: Thermal voltage for the module [V]
        - Rs, Rsh: Resistances [Ω]
    """
    T_kelvin = T_celsius + 273.15
    T_stc_kelvin = module.T_stc + 273.15
    dT = T_celsius - module.T_stc
    
    # Thermal voltage for entire module (Ns_cells in series)
    Vt = module.n * module.Ns_cells * K_BOLTZMANN * T_kelvin / Q_ELECTRON
    Vt_stc = module.n * module.Ns_cells * K_BOLTZMANN * T_stc_kelvin / Q_ELECTRON
    
    # Light-generated current scales linearly with irradiance and with temperature
    IL = (G / module.G_stc) * (module.Isc_stc + module.alpha_Isc * dT)
    
    # Saturation current: Io(T) = Io(T_stc) * (T/T_stc)^3 * exp(qEg/nk * (1/T_stc - 1/T))
    # Simplified: we back-calculate Io from STC conditions
    # At STC, near open circuit: Voc ≈ Vt * ln(IL/Io)
    Io_stc = module.Isc_stc / (np.exp(module.Voc_stc / Vt_stc) - 1)
    
    # Temperature scaling of Io (approximate, using exponential factor)
    # Using Eg ≈ 1.12 eV for silicon
    Eg = 1.12  # eV
    Io = Io_stc * (T_kelvin / T_stc_kelvin)**3 * np.exp(
        (Q_ELECTRON * Eg / (module.n * K_BOLTZMANN)) * (1/T_stc_kelvin - 1/T_kelvin)
    )
    
    return {
        'IL': IL,
        'Io': Io,
        'Vt': Vt,
        'Rs': module.Rs,
        'Rsh': module.Rsh
    }


def _substring_voltage_at_current(
    I: float,
    IL: float,
    Io: float,
    Rs: float,
    Rsh: float,
    Vt: float,
    max_iter: int = 50,
    tol: float = 1e-6
) -> float:
    """
    Compute substring voltage for a given current using Newton-Raphson.
    
    Single-diode equation for a substring:
    I = IL - Io * (exp((V + I*Rs)/Vt) - 1) - (V + I*Rs)/Rsh
    
    Solve for V given I using Newton-Raphson on f(V) = 0:
    f(V) = I - IL + Io * (exp((V + I*Rs)/Vt) - 1) + (V + I*Rs)/Rsh
    
    Args:
        I: Current [A]
        IL, Io, Rs, Rsh, Vt: Scaled single-diode parameters for substring
    
    Returns:
        Substring voltage [V]
    """
    # Initial guess: V ≈ -I*Rs (works well for most cases)
    V = -I * Rs
    
    for iteration in range(max_iter):
        # Compute exponent with clamping for numerical stability
        exponent = (V + I * Rs) / Vt
        exponent = np.clip(exponent, -40, 40)  # Prevent overflow
        exp_term = np.exp(exponent)
        
        # Function value
        f = I - IL + Io * (exp_term - 1) + (V + I * Rs) / Rsh
        
        # Derivative df/dV
        df_dV = Io * exp_term / Vt + 1 / Rsh
        
        # Guard against zero derivative
        if abs(df_dV) < 1e-12:
            break
        
        # Newton update
        dV = -f / df_dV
        V = V + dV
        
        # Check convergence
        if abs(dV) < tol and abs(f) < tol:
            break
    
    return V


def module_voltage_at_current_with_substrings(
    I: float,
    G: float,
    T_celsius: float,
    module: ModuleParams
) -> float:
    """
    Compute module voltage at a given current, accounting for bypass diodes.
    
    The module has 3 substrings, each with 18 cells. For the given series
    current I, we compute each substring's voltage. If a substring would
    operate below -Vd (bypass diode forward drop), we clamp it to
    V_substring = -(Vd + I * rd) to model bypass conduction.
    
    Args:
        I: Module current [A]
        G: Irradiance on this module [W/m²]
        T_celsius: Cell temperature [°C]
        module: Module parameters
    
    Returns:
        Module voltage [V] (sum of 3 substring voltages)
    """
    if G < 1e-3:  # Very low irradiance, module is essentially off
        return -3 * (module.Vd + I * module.rd) if I > 0 else 0.0
    
    # Scale parameters for operating conditions
    params = scale_params(module, G, T_celsius)
    
    # Scale to per-substring parameters (18 cells out of 54 total)
    substring_fraction = module.Ns_per_substring / module.Ns_cells
    IL_sub = params['IL']  # Current is same
    Io_sub = params['Io']  # Saturation current is same
    Rs_sub = params['Rs'] * substring_fraction  # Resistance scales
    Rsh_sub = params['Rsh'] * substring_fraction
    Vt_sub = params['Vt'] * substring_fraction  # Thermal voltage scales
    
    # Compute voltage for each substring
    V_module = 0.0
    bypass_threshold = -module.Vd
    
    for substring_idx in range(3):
        # Compute ideal substring voltage
        V_sub = _substring_voltage_at_current(
            I, IL_sub, Io_sub, Rs_sub, Rsh_sub, Vt_sub
        )
        
        # Check if bypass diode conducts
        if V_sub < bypass_threshold:
            # Bypass diode is conducting, clamp voltage
            V_sub = -(module.Vd + I * module.rd)
        
        V_module += V_sub
    
    return V_module


def solve_parallel_block_voltage_for_current(
    I_target: float,
    G_row: np.ndarray,
    T_row: np.ndarray,
    module: ModuleParams,
    max_iter: int = 30,
    tol: float = 1e-4
) -> float:
    """
    Find the voltage Vr across a parallel block (row) for a target current.
    
    In TCT topology, each row is a parallel connection of C modules. All modules
    in the row share the same voltage Vr, but carry different currents depending
    on their irradiance. The sum of module currents must equal I_target.
    
    This inverts the module V(I) relationship to find Vr such that:
        sum_over_modules( I_module(Vr, G_c, T_c) ) = I_target
    
    We allow small negative Vr (down to ~-1 V) to model bypass conduction.
    
    Args:
        I_target: Target row current [A]
        G_row: Irradiance for each module in row [W/m²], shape (C,)
        T_row: Temperature for each module in row [°C], shape (C,)
        module: Module parameters
    
    Returns:
        Row voltage Vr [V]
    """
    C = len(G_row)  # Number of modules in parallel
    
    # Initial guess: use typical operating voltage
    Vr = module.Vmpp_stc
    
    # Helper function to compute total row current at voltage Vr
    def row_current_at_voltage(Vr: float) -> Tuple[float, float]:
        """Returns (I_total, dI_total/dVr)"""
        I_total = 0.0
        dI_dVr = 0.0
        
        for c in range(C):
            G_c = G_row[c]
            T_c = T_row[c]
            
            if G_c < 1e-3:  # Module is essentially off
                continue
            
            # Scale parameters
            params = scale_params(module, G_c, T_c)
            IL = params['IL']
            Io = params['Io']
            Vt = params['Vt']
            Rs = params['Rs']
            Rsh = params['Rsh']
            
            # Single-diode equation: I = IL - Io*(exp((V+I*Rs)/Vt)-1) - (V+I*Rs)/Rsh
            # Solve for I given V=Vr using Newton-Raphson
            I_c = IL * 0.9  # Initial guess
            
            for _ in range(20):
                exponent = np.clip((Vr + I_c * Rs) / Vt, -40, 40)
                exp_term = np.exp(exponent)
                
                f = I_c - IL + Io * (exp_term - 1) + (Vr + I_c * Rs) / Rsh
                df_dI = 1 + Io * Rs * exp_term / Vt + Rs / Rsh
                
                if abs(df_dI) < 1e-12:
                    break
                
                dI = -f / df_dI
                I_c = I_c + dI
                
                if abs(dI) < 1e-6:
                    break
            
            # Derivative dI/dVr (implicit differentiation)
            exponent = np.clip((Vr + I_c * Rs) / Vt, -40, 40)
            exp_term = np.exp(exponent)
            df_dV = Io * exp_term / Vt + 1 / Rsh
            df_dI = 1 + Io * Rs * exp_term / Vt + Rs / Rsh
            dI_c_dVr = -df_dV / df_dI if abs(df_dI) > 1e-12 else 0
            
            I_total += I_c
            dI_dVr += dI_c_dVr
        
        return I_total, dI_dVr
    
    # Newton-Raphson to find Vr such that row_current(Vr) = I_target
    for iteration in range(max_iter):
        I_current, dI_dVr = row_current_at_voltage(Vr)
        
        error = I_current - I_target
        
        # Guard against zero derivative
        if abs(dI_dVr) < 1e-12:
            break
        
        # Newton update
        dVr = -error / dI_dVr
        
        # Limit step size
        dVr = np.clip(dVr, -5, 5)
        
        Vr = Vr + dVr
        
        # Allow small negative voltage for bypass conduction
        Vr = max(Vr, -1.0)
        
        if abs(error) < tol and abs(dVr) < tol:
            break
    
    return Vr


def evaluate_tct(
    G_map: np.ndarray,
    T_map: np.ndarray,
    module: ModuleParams,
    I_points: Optional[np.ndarray] = None,
    num_points: int = 500
) -> Dict:
    """
    Evaluate I-V and P-V curves for a TCT-connected PV array.
    
    TCT topology:
    - Each row is a parallel block of C modules (all at same voltage)
    - Rows are connected in series (same current through all rows)
    
    Algorithm:
    1. Find array Isc by solving V_array(I) = 0
    2. Sweep current from 0 to Isc
    3. For each current I, find each row's voltage Vr(I)
    4. Array voltage = sum of row voltages
    5. Find GMPP and local MPPs
    
    Args:
        G_map: Irradiance map [W/m²], shape (R, C) where R=rows, C=columns
        T_map: Temperature map [°C], shape (R, C)
        module: Module parameters
        I_points: Optional array of current points to evaluate [A]
        num_points: Number of points if I_points not provided
    
    Returns:
        Dictionary with:
        - 'I': Current array [A]
        - 'V': Voltage array [V]
        - 'P': Power array [W]
        - 'gmpp': Tuple (Pmpp, Impp, Vmpp) for global MPP
        - 'local_mpps': List of (P, I, V) tuples for all local maxima
    """
    R, C = G_map.shape
    
    # Find array short-circuit current (V=0)
    # For V=0, each row voltage is 0, so we need to find max current
    def array_voltage_at_current(I: float) -> float:
        """Compute total array voltage for given current"""
        V_total = 0.0
        for r in range(R):
            Vr = solve_parallel_block_voltage_for_current(
                I, G_map[r], T_map[r], module
            )
            V_total += Vr
        return V_total
    
    # Find Isc using bisection on V(I) = 0
    # Upper bound: sum of all module Isc values (conservative)
    I_max = np.sum(G_map) / module.G_stc * module.Isc_stc
    
    # Bisection search
    I_low, I_high = 0.0, I_max
    for _ in range(50):
        I_mid = (I_low + I_high) / 2
        V_mid = array_voltage_at_current(I_mid)
        
        if V_mid > 0:
            I_low = I_mid
        else:
            I_high = I_mid
        
        if I_high - I_low < 1e-4:
            break
    
    Isc_array = I_high
    
    # Generate current sweep points
    if I_points is None:
        # Dense sampling near Isc for better resolution
        I_points = np.concatenate([
            np.linspace(0, Isc_array * 0.95, num_points - 50),
            np.linspace(Isc_array * 0.95, Isc_array, 50)
        ])
    
    # Evaluate array voltage and power at each current point
    V_array = np.zeros_like(I_points)
    
    for i, I in enumerate(I_points):
        V_array[i] = array_voltage_at_current(I)
    
    P_array = V_array * I_points
    
    # Find global MPP
    idx_gmpp = np.argmax(P_array)
    gmpp = (P_array[idx_gmpp], I_points[idx_gmpp], V_array[idx_gmpp])
    
    # Find local MPPs
    local_mpps = find_local_mpps(I_points, V_array, P_array)
    
    return {
        'I': I_points,
        'V': V_array,
        'P': P_array,
        'gmpp': gmpp,
        'local_mpps': local_mpps
    }


def find_local_mpps(
    I: np.ndarray,
    V: np.ndarray,
    P: np.ndarray,
    prominence_threshold: float = 0.01
) -> List[Tuple[float, float, float]]:
    """
    Find local maximum power points in the P-V curve.
    
    A local MPP is identified where dP/dI changes from positive to negative,
    and the peak is sufficiently prominent (not just noise).
    
    Args:
        I: Current array [A]
        V: Voltage array [V]
        P: Power array [W]
        prominence_threshold: Minimum relative prominence (fraction of GMPP)
    
    Returns:
        List of (P, I, V) tuples for local maxima, sorted by power (descending)
    """
    local_mpps = []
    
    if len(P) < 3:
        return local_mpps
    
    # Compute derivative dP/dI using central differences
    dP_dI = np.gradient(P, I)
    
    # Find zero-crossings from positive to negative
    for i in range(1, len(dP_dI) - 1):
        if dP_dI[i-1] > 0 and dP_dI[i+1] < 0:
            # Local maximum at index i
            P_local = P[i]
            
            # Check prominence (compare to GMPP)
            P_gmpp = np.max(P)
            if P_local > prominence_threshold * P_gmpp:
                local_mpps.append((P_local, I[i], V[i]))
    
    # Sort by power (descending)
    local_mpps.sort(reverse=True)
    
    return local_mpps


def generate_demo_shading(
    R: int,
    C: int,
    levels: List[int] = [100, 200, 400, 600, 800, 1000],
    pattern: str = 'diagonal'
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate demonstration shading patterns for testing.
    
    Args:
        R: Number of rows
        C: Number of columns
        levels: Irradiance levels [W/m²]
        pattern: Shading pattern type ('diagonal', 'checkerboard', 'random', etc.)
    
    Returns:
        Tuple of (G_map, T_map):
        - G_map: Irradiance map [W/m²], shape (R, C)
        - T_map: Temperature map [°C], shape (R, C) - all at 25°C for simplicity
    """
    G_map = np.zeros((R, C))
    
    if pattern == 'diagonal':
        # Diagonal gradient pattern
        for r in range(R):
            for c in range(C):
                idx = min((r + c) * len(levels) // (R + C), len(levels) - 1)
                G_map[r, c] = levels[idx]
    
    elif pattern == 'checkerboard':
        # Checkerboard pattern
        for r in range(R):
            for c in range(C):
                idx = ((r % 2) + (c % 2)) % 2
                G_map[r, c] = levels[idx * (len(levels) - 1)]
    
    elif pattern == 'random':
        # Random pattern
        G_map = np.random.choice(levels, size=(R, C))
    
    elif pattern == 'partial':
        # Partial shading: top-left corner shaded
        G_map.fill(1000)
        shade_R = R // 2
        shade_C = C // 2
        for r in range(shade_R):
            for c in range(shade_C):
                G_map[r, c] = levels[0]
    
    else:
        # Uniform (no shading)
        G_map.fill(1000)
    
    # Temperature uniform at 25°C
    T_map = np.full((R, C), 25.0)
    
    return G_map, T_map


def plot_iv_pv(
    results: Dict,
    title: str = "TCT Array I-V and P-V Curves",
    show: bool = True,
    filename: Optional[str] = None
):
    """
    Plot I-V and P-V curves with MPP markers.
    
    Args:
        results: Dictionary from evaluate_tct()
        title: Plot title
        show: Whether to display the plot
        filename: Optional filename to save the plot
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Matplotlib not available. Install it to use plotting: pip install matplotlib")
        return
    
    I = results['I']
    V = results['V']
    P = results['P']
    gmpp = results['gmpp']
    local_mpps = results['local_mpps']
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # I-V curve
    ax1.plot(V, I, 'b-', linewidth=2, label='I-V Curve')
    ax1.plot(gmpp[2], gmpp[1], 'r*', markersize=15, label=f'GMPP ({gmpp[0]:.1f} W)')
    
    for i, (P_local, I_local, V_local) in enumerate(local_mpps[1:]):
        ax1.plot(V_local, I_local, 'go', markersize=8, 
                label=f'Local MPP {i+1} ({P_local:.1f} W)' if i < 3 else '')
    
    ax1.set_xlabel('Voltage [V]', fontsize=12)
    ax1.set_ylabel('Current [A]', fontsize=12)
    ax1.set_title('I-V Characteristic', fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # P-V curve
    ax2.plot(V, P, 'b-', linewidth=2, label='P-V Curve')
    ax2.plot(gmpp[2], gmpp[0], 'r*', markersize=15, label=f'GMPP ({gmpp[0]:.1f} W)')
    
    for i, (P_local, I_local, V_local) in enumerate(local_mpps[1:]):
        ax2.plot(V_local, P_local, 'go', markersize=8,
                label=f'Local MPP {i+1} ({P_local:.1f} W)' if i < 3 else '')
    
    ax2.set_xlabel('Voltage [V]', fontsize=12)
    ax2.set_ylabel('Power [W]', fontsize=12)
    ax2.set_title('P-V Characteristic', fontsize=14)
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    fig.suptitle(title, fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    if filename:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {filename}")
    
    if show:
        plt.show()
    else:
        plt.close()


# Example usage and testing
if __name__ == "__main__":
    print("=" * 70)
    print("TCT PV Array Evaluation Module - Demo")
    print("=" * 70)
    
    # Create module
    module = default_kc200gt()
    print(f"\nModule: KC200GT")
    print(f"  Voc = {module.Voc_stc} V, Isc = {module.Isc_stc} A")
    print(f"  Vmpp = {module.Vmpp_stc} V, Impp = {module.Impp_stc} A")
    print(f"  Pmpp = {module.Vmpp_stc * module.Impp_stc:.1f} W")
    
    # Test 1: Uniform irradiance (no shading)
    print("\n" + "-" * 70)
    print("Test 1: Uniform irradiance (4x3 array, 1000 W/m²)")
    print("-" * 70)
    
    R, C = 4, 3
    G_uniform, T_uniform = generate_demo_shading(R, C, pattern='uniform')
    
    results_uniform = evaluate_tct(G_uniform, T_uniform, module, num_points=300)
    
    gmpp = results_uniform['gmpp']
    print(f"\nGlobal MPP:")
    print(f"  Pmpp = {gmpp[0]:.2f} W")
    print(f"  Impp = {gmpp[1]:.2f} A")
    print(f"  Vmpp = {gmpp[2]:.2f} V")
    print(f"  Isc = {results_uniform['I'][-1]:.2f} A")
    print(f"  Voc = {results_uniform['V'][0]:.2f} V")
    
    # Test 2: Partial shading
    print("\n" + "-" * 70)
    print("Test 2: Partial shading (4x3 array, corner shaded)")
    print("-" * 70)
    
    G_partial, T_partial = generate_demo_shading(R, C, pattern='partial')
    print(f"\nIrradiance map [W/m²]:")
    print(G_partial)
    
    results_partial = evaluate_tct(G_partial, T_partial, module, num_points=300)
    
    gmpp = results_partial['gmpp']
    print(f"\nGlobal MPP:")
    print(f"  Pmpp = {gmpp[0]:.2f} W")
    print(f"  Impp = {gmpp[1]:.2f} A")
    print(f"  Vmpp = {gmpp[2]:.2f} V")
    
    print(f"\nLocal MPPs found: {len(results_partial['local_mpps'])}")
    for i, (P, I, V) in enumerate(results_partial['local_mpps'][:5]):
        print(f"  MPP {i+1}: P={P:.2f} W, I={I:.2f} A, V={V:.2f} V")
    
    # Test 3: Diagonal shading
    print("\n" + "-" * 70)
    print("Test 3: Diagonal shading pattern (4x3 array)")
    print("-" * 70)
    
    G_diagonal, T_diagonal = generate_demo_shading(
        R, C, pattern='diagonal', levels=[200, 400, 600, 800, 1000]
    )
    print(f"\nIrradiance map [W/m²]:")
    print(G_diagonal)
    
    results_diagonal = evaluate_tct(G_diagonal, T_diagonal, module, num_points=300)
    
    gmpp = results_diagonal['gmpp']
    print(f"\nGlobal MPP:")
    print(f"  Pmpp = {gmpp[0]:.2f} W")
    print(f"  Impp = {gmpp[1]:.2f} A")
    print(f"  Vmpp = {gmpp[2]:.2f} V")
    
    print(f"\nLocal MPPs found: {len(results_diagonal['local_mpps'])}")
    for i, (P, I, V) in enumerate(results_diagonal['local_mpps'][:5]):
        print(f"  MPP {i+1}: P={P:.2f} W, I={I:.2f} A, V={V:.2f} V")
    
    # Generate plots if matplotlib is available
    print("\n" + "-" * 70)
    print("Generating plots...")
    print("-" * 70)

    try:
        plot_iv_pv(results_uniform, title="Uniform Irradiance (No Shading)", 
                show=False, filename="tct_uniform.png")
        plot_iv_pv(results_partial, title="Partial Shading", 
                show=False, filename="tct_partial.png")
        plot_iv_pv(results_diagonal, title="Diagonal Shading Pattern", 
                show=False, filename="tct_diagonal.png")
        print("Plots saved successfully!")

    except Exception as e:
        print(f"Could not generate plots: {e}")
    
    print("\n" + "=" * 70)
    print("Demo completed successfully!")
    print("=" * 70)
