# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Run the built-in demo (uniform, partial, and diagonal shading scenarios)
python tct_eval.py

# Run the verification test suite (8 scenarios, validates numerical accuracy)
python verify_tct.py

# Run array-specific test scripts
python test_3x4_array.py
python test_9x9_array.py
```

Dependencies: `numpy` (required), `matplotlib` (optional, for plots).

## Architecture

This repository is a single-module Python library (`tct_eval.py`) for simulating PV arrays under partial shading. There is no build system, package structure, or test framework — scripts import directly from `tct_eval.py` in the same directory.

### Core abstraction: `ModuleParams` dataclass

All physical parameters for a PV module live in `ModuleParams`. `default_kc200gt()` returns a pre-populated instance for the Kyocera KC200GT (200W, 54 cells, 3 substrings × 18 cells).

### Computation pipeline in `evaluate_tct(G_map, T_map, module)`

The main entry point. `G_map` and `T_map` are `(R, C)` numpy arrays where R = rows (series) and C = columns (parallel).

1. **Isc search** — bisection on `array_voltage_at_current(I) = 0`
2. **Current sweep** — for each current point I:
   - Per row: `solve_parallel_block_voltage_for_current()` — Newton-Raphson to find row voltage Vr such that sum of module currents = I
   - Per module in row: `module_voltage_at_current_with_substrings()` — solves the single-diode equation for each of 3 substrings via `_substring_voltage_at_current()` (Newton-Raphson), then clamps any substring below −Vd to model bypass diode conduction
3. **MPP detection** — `find_local_mpps()` uses `np.gradient` on P(I) and finds zero-crossings of dP/dI

### TCT topology convention

- **Rows** (axis 0) are connected in **series** — same current, voltages add
- **Columns** (axis 1) are connected in **parallel** — same voltage, currents add
- `G_map[r, c]` = irradiance on the module at row r, column c

### Results dictionary

`evaluate_tct` returns `{'I', 'V', 'P', 'gmpp', 'local_mpps'}` where `gmpp` is `(Pmpp, Impp, Vmpp)` and `local_mpps` is a list of the same tuple sorted by power descending.

### Numerical notes

- Exponential terms are clamped to `[-40, 40]` to prevent overflow
- Newton-Raphson steps are bounded to ±5 V at the row solver level
- Row voltage is floor-clamped at −1 V to allow bypass conduction without diverging
- The `np.gradient` divide-by-zero RuntimeWarnings at startup are benign — they occur at array boundaries during MPP detection and do not affect results
