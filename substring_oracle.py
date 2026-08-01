"""Substring-level TCT cross-validation oracle.

Promoted to the regression-suite oracle in the round-3 Phase 2a fix
(2026-07-30, decision R1) after the Phase 1 validation gate (T3a) showed
the -1.0 V row clamp shared by the pre-fix backend and the monolithic
``tct_eval.evaluate_tct`` path relocates the global maximum power point
on one released pattern. ``tct_eval.py`` is retained beside this file
for provenance and for its single-diode primitives, which this oracle
reuses. Original gate copy archived under the manuscript folder at
``code/2026-07-30-phase1-pv-gate/substring_tct.py``.

Substring-level TCT solver, built for the round-3 Phase 1 PV validation gate (T3a).

Independent comparison model built from the cross-validation oracle's unused
per-substring machinery (tct_eval.py: scale_params,
_substring_voltage_at_current, module_voltage_at_current_with_substrings,
round-2 finding B17).

Physical model
--------------
Each KC200GT module is composed of three 18-cell substrings. Each substring is
a single-diode device with parameters partitioned from the module parameters
(Rs, Rsh, and thermal voltage scaled by 18/54), and each substring is
protected by its own bypass diode modelled as a 0.5 V forward drop in series
with a 0.02 Ohm dynamic resistance (registry values). A substring whose
single-diode voltage at the imposed current would fall below -Vd is clamped to
V = -(Vd + I*rd), i.e. its bypass diode conducts. Module voltage is the sum of
the three substring voltages. Rows are parallel blocks of substring-modelled
modules (shared row voltage, currents sum), and rows are connected in series
(shared current, voltages sum). There is NO row-level -1.0 V clamp: row
voltage is bounded only by the physical bypass composition.

Constraint stated per the register: every module receives spatially uniform
irradiance internally (one irradiance value per module), so all three
substrings of a module share the same G and T. Under uniform G within the
module, the three substrings are identical; the substring model therefore
differs from the released monolithic-module model only where bypass
conduction is engaged (reverse-biased modules inside a row, reverse-biased
rows), because in the forward region the series composition of three 18-cell
substrings is algebraically identical to one 54-cell device.

Numerics
--------
Per unique (G, T) pair, the substring voltage V_sub(I) is computed on a dense
current grid by vectorised Newton-Raphson (same equation and parameter
scaling as the oracle; verified pointwise against the oracle's scalar
module_voltage_at_current_with_substrings). The module curve is
V_mod(I) = sum over substrings of max(V_sub(I), -(Vd + I*rd)).
The current grid is globally dense and additionally refined around each
unique photocurrent value, where the bypass knee is sharp. Module curves are
monotone decreasing in I, so I_mod(V) is obtained by inverse interpolation.
Row current at voltage Vr is the sum of module currents; row voltage at a
target current is found by bisection on the monotone row current function.
Array voltage at a current is the sum of row voltages; array Isc is found by
bisection on V_array(I) = 0, and the I sweep then mirrors the released
solver's construction (num_points with a dense tail above 0.95*Isc).
"""

import sys
from pathlib import Path

import numpy as np

# tct_eval lives in the same directory; make the import work whether this
# file is imported as a module or run as a script.
_HERE = str(Path(__file__).resolve().parent)
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import tct_eval
from tct_eval import ModuleParams, default_kc200gt, scale_params


def substring_voltage_vec(I_grid, G, T, module, max_iter=300, tol=1e-12):
    """Vectorised Newton-Raphson for substring voltage at each current.

    Same single-diode equation and parameter partitioning as the oracle's
    _substring_voltage_at_current / module_voltage_at_current_with_substrings.
    Two numerical robustness additions over the oracle scalar (which is only
    exercised at I in [0, ~Isc]): a diode-informed initial guess valid for
    reverse currents (I < 0, V > Voc) and deep forward overload (I >> Isc),
    and a +/-2 V Newton step clamp to prevent overshoot through the clipped
    exponential. Verified pointwise against the oracle scalar in
    verify_against_oracle_scalar().
    """
    params = scale_params(module, G, T)
    frac = module.Ns_per_substring / module.Ns_cells
    IL = params["IL"]
    Io = params["Io"]
    Rs = params["Rs"] * frac
    Rsh = params["Rsh"] * frac
    Vt = params["Vt"] * frac
    # Initial guess: diode branch where I < IL (V ~ Vt*ln((IL-I)/Io+1) - I*Rs),
    # shunt branch where I >= IL (V ~ Rsh*(IL-I) - I*Rs).
    diode_arg = np.maximum((IL - I_grid) / Io, 0.0) + 1.0
    V_diode = Vt * np.log(diode_arg) - I_grid * Rs
    V_shunt = Rsh * (IL - I_grid) - I_grid * Rs
    V = np.where(I_grid < IL, V_diode, V_shunt)
    for _ in range(max_iter):
        expo = np.clip((V + I_grid * Rs) / Vt, -40.0, 40.0)
        e = np.exp(expo)
        f = I_grid - IL + Io * (e - 1.0) + (V + I_grid * Rs) / Rsh
        df = Io * e / Vt + 1.0 / Rsh
        dV = np.clip(-f / df, -2.0, 2.0)
        V = V + dV
        if np.max(np.abs(dV)) < tol:
            break
    return V


def module_voltage_curve(I_grid, G, T, module):
    """Module V(I) with per-substring bypass diodes (three substrings).

    Uniform irradiance within the module: the three substrings share G and
    T, so one substring solve serves all three. Bypass clamp applied per
    substring at V = -(Vd + I*rd).
    """
    if G < 1e-3:
        # Oracle convention for a dark module.
        return np.where(I_grid > 0, -3.0 * (module.Vd + I_grid * module.rd), 0.0)
    V_sub = substring_voltage_vec(I_grid, G, T, module)
    clamp = -(module.Vd + I_grid * module.rd)
    return 3.0 * np.maximum(V_sub, clamp)


def _build_I_grid(G_map, module, n_base=40001, n_knee=4001):
    """Dense current grid with refinement around each bypass knee."""
    IL_vals = sorted(set((G / module.G_stc) * module.Isc_stc
                         for G in np.asarray(G_map).ravel() if G >= 1e-3))
    row_sums = [np.sum(row) / module.G_stc * module.Isc_stc
                for row in np.asarray(G_map)]
    I_hi = 1.05 * max(row_sums)
    I_lo = -0.35 * module.Isc_stc  # allow reverse module current near row Voc
    pieces = [np.linspace(I_lo, I_hi, n_base)]
    for IL in IL_vals:
        pieces.append(np.linspace(max(I_lo, IL - 0.08), min(I_hi, IL + 0.08), n_knee))
    grid = np.unique(np.concatenate(pieces))
    return grid


class SubstringTCTModel:
    """Substring-level TCT array model (rows parallel, rows series)."""

    def __init__(self, G_map, T_map, module=None):
        self.module = module if module is not None else default_kc200gt()
        self.G = np.asarray(G_map, dtype=float)
        self.T = np.asarray(T_map, dtype=float)
        self.R, self.C = self.G.shape
        self.I_grid = _build_I_grid(self.G, self.module)
        # Per unique (G, T): module curve on the grid.
        self._curves = {}
        for r in range(self.R):
            for c in range(self.C):
                key = (float(self.G[r, c]), float(self.T[r, c]))
                if key not in self._curves:
                    self._curves[key] = module_voltage_curve(
                        self.I_grid, key[0], key[1], self.module)
        # Per-module inverted curves (V ascending -> I) for row summation.
        self._inv = {}
        for key, Vc in self._curves.items():
            # V(I) monotone decreasing: reverse for ascending interp input.
            self._inv[key] = (Vc[::-1].copy(), self.I_grid[::-1].copy())
        # Row voltage search bounds.
        self._Vr_lo = float(min(Vc.min() for Vc in self._curves.values()))
        self._Vr_hi = float(max(Vc.max() for Vc in self._curves.values()))

    def row_current(self, Vr, r):
        """Sum of module currents in row r at row voltage Vr."""
        total = 0.0
        for c in range(self.C):
            key = (float(self.G[r, c]), float(self.T[r, c]))
            V_asc, I_desc = self._inv[key]
            total += float(np.interp(Vr, V_asc, I_desc))
        return total

    def row_voltage_at_current(self, I_target, r, iters=80):
        """Bisection on the monotone-decreasing row current function."""
        lo, hi = self._Vr_lo, self._Vr_hi
        # row_current is decreasing in Vr: current(lo) is max, current(hi) is min.
        if I_target >= self.row_current(lo, r):
            return lo
        if I_target <= self.row_current(hi, r):
            return hi
        for _ in range(iters):
            mid = 0.5 * (lo + hi)
            if self.row_current(mid, r) > I_target:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    def array_voltage_at_current(self, I_target):
        return sum(self.row_voltage_at_current(I_target, r)
                   for r in range(self.R))

    def sweep(self, num_points=500):
        """I-V/P-V sweep mirroring the released solver's construction."""
        # Bisection for array Isc on V_array(I) = 0 (released convention).
        I_max = float(self.G.sum() / self.module.G_stc * self.module.Isc_stc)
        I_lo_b, I_hi_b = 0.0, I_max
        for _ in range(60):
            I_mid = 0.5 * (I_lo_b + I_hi_b)
            if self.array_voltage_at_current(I_mid) > 0:
                I_lo_b = I_mid
            else:
                I_hi_b = I_mid
            if I_hi_b - I_lo_b < 1e-6:
                break
        Isc = I_hi_b
        n_dense = min(50, max(num_points // 10, 10))
        n_coarse = num_points - n_dense
        I_pts = np.concatenate([
            np.linspace(0.0, Isc * 0.95, n_coarse),
            np.linspace(Isc * 0.95, Isc, n_dense),
        ])
        V_pts = np.array([self.array_voltage_at_current(float(I))
                          for I in I_pts])
        P_pts = V_pts * I_pts
        idx = int(np.argmax(P_pts))
        return {
            "I": I_pts, "V": V_pts, "P": P_pts,
            "Pmpp": float(P_pts[idx]), "Vmpp": float(V_pts[idx]),
            "Impp": float(I_pts[idx]), "Voc": float(V_pts[0]),
            "Isc": float(Isc),
        }


def verify_solution(model, n_samples=40, seed=0):
    """Three-layer verification of the module curves.

    1. Residual check: the single-diode residual f(V_sub) at every grid
       point in the non-bypassed region must be ~0 (the clamped region is
       exact by construction).
    2. Independent root-find: scipy.optimize.brentq on the same residual
       (f is strictly increasing in V, so the root is unique) at random
       currents, compared to grid interpolation.
    3. Oracle-scalar agreement near the bypass knee, the region where the
       oracle's unused module_voltage_at_current_with_substrings converges.
       (Outside that region the oracle scalar diverges: its Newton starts
       at V = -I*Rs and takes an unclamped first step of |f|/f' ~ hundreds
       of volts into the clipped exponential, then cannot walk back within
       its 50-iteration budget. This latent defect is why layers 1-2, not
       the oracle scalar, are the primary verification.)

    Returns dict of max abs errors per layer.
    """
    from scipy.optimize import brentq
    rng = np.random.default_rng(seed)
    p = model.module
    frac = p.Ns_per_substring / p.Ns_cells
    out = {"residual": 0.0, "brentq": 0.0, "oracle_knee": 0.0}
    for (G, T), Vmod in model._curves.items():
        if G < 1e-3:
            continue
        sp = scale_params(p, G, T)
        IL, Io = sp["IL"], sp["Io"]
        Rs, Rsh, Vt = sp["Rs"] * frac, sp["Rsh"] * frac, sp["Vt"] * frac

        def f(V, I):
            return I - IL + Io * (np.exp(np.clip((V + I * Rs) / Vt, -40, 40)) - 1.0) + (V + I * Rs) / Rsh

        # Layer 1: residual at grid points not clamped by the bypass model.
        V_sub = Vmod / 3.0
        clamp = -(p.Vd + model.I_grid * p.rd)
        free = V_sub > clamp + 1e-9
        res = np.abs(f(V_sub[free], model.I_grid[free]))
        out["residual"] = max(out["residual"], float(res.max()))
        # Layer 2: brentq at random currents across the full sweep range.
        I_samples = rng.uniform(model.I_grid.min(), model.I_grid.max(), n_samples)
        for I in I_samples:
            lo = min(-2.0, Rsh * (IL - float(I)) - float(I) * Rs - 5.0)
            v_root = brentq(f, lo, 45.0, args=(float(I),), xtol=1e-12)
            v_model = max(v_root, -(p.Vd + I * p.rd))
            v_grid = float(np.interp(I, model.I_grid, Vmod)) / 3.0
            out["brentq"] = max(out["brentq"], abs(v_model - v_grid))
        # Layer 3: oracle scalar near the knee (its convergent region).
        for dI in (-0.02, 0.0, 0.02, 0.2, 0.5):
            I = IL + dI
            v_oracle = tct_eval.module_voltage_at_current_with_substrings(
                float(I), G, T, p)
            v_grid = float(np.interp(I, model.I_grid, Vmod))
            out["oracle_knee"] = max(out["oracle_knee"], abs(v_oracle - v_grid))
    return out
