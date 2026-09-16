# zBET

**zBET** is an advanced rotor aerodynamic solver based on Blade Element Momentum Theory (BEMT). It computes rotor trim, thrust, shaft torque, longitudinal forces, rolling and pitching moments, and power metrics across advance ratio sweeps ($\mu$) for helicopters and rotary-wing aircraft.

For the complete theoretical derivations, mathematical formulas, and literature citations, see the [zBET Documentation](zBET-documentation.md).

---

## Features

- **4 Inflow Models**:
  - Classic Uniform Inflow (Glauert, 1926)
  - Simple Coleman Inflow (Coleman et al., 1945)
  - NDARC Coleman-Feingold Inflow with lateral gradient $K_y$ (NASA/TP-2009-215402)
  - Drees Inflow with longitudinal and lateral gradients (Drees, 1949)
- **3 Solidity Options & Taper**:
  - Reference solidity $\sigma_{\mathrm{ref}}$ (extrapolated to hub $r=0$)
  - Geometric physical solidity $\sigma_{\mathrm{geom}}$ (actual blade area / disk area)
  - Thrust-weighted solidity $\sigma_{\mathrm{thrust}}$ ($r^2$-weighted)
  - Linear blade chord taper ($c_{\mathrm{root}}$ to $c_{\mathrm{tip}}$)
- **Blade Pitch & Hover Trim Modes**:
  - Collective trim with invariant total twist ($\Delta\theta = \text{const}$)
  - RPM trim targeting rotor thrust in Newtons ($T$)
  - Constant pitch or linear twist ($\theta_{\mathrm{root}}$ to $\theta_{\mathrm{tip}}$)
- **Advanced Aerodynamic Corrections**:
  - Blade tip loss factor $B$ (fixed or Sissingh self-adjusting formula)
  - Prandtl-Glauert compressibility correction $a(M)$ on section lift curve slope
- **Comprehensive Outputs**:
  - Non-dimensional coefficients: $C_T, C_Q, C_{Qi}, C_{Q0}, C_H, C_{Hi}, C_{H0}, C_Y, C_{Mx}, C_{My}$
  - Performance metrics: Total air power $C_{P,\mathrm{air}} = C_Q + \mu C_H$, effective rotor lift-to-drag $(L/D)_{\mathrm{eff}}$, and hover Figure of Merit ($FoM$)
  - Dimensional outputs: Thrust $T$ [N], Shaft Power $P$ [kW], and Total Air Power $P_{\mathrm{air}}$ [kW]

---

## Quick Start

### Requirements
- Python 3.9+
- `numpy`, `pandas`, `matplotlib`, `pytest`

Install dependencies:
```bash
pip install numpy pandas matplotlib pytest
```

### Running the Solver
Run the main sweep using standard utility helicopter parameters:
```bash
python zBET.py
```

Outputs are automatically generated in the `outputs/` directory:
- `outputs/zBET.csv` (consolidated dataset with all inflow models)
- `outputs/zBET_<model>.csv` (individual CSV per inflow model)
- 14 performance plots vs. $\mu$ for each inflow model

### Running Tests
Run the automated test suite:
```bash
pytest
```

---

## Project Structure
```text
zBET/
├── zBET.py                  # Core BEMT solver and sweep engine
├── zBET-documentation.md    # In-depth mathematical theory and derivations
├── README.md                # Project overview and quickstart guide
├── AGENTS.md                # AI agent developer guidelines
├── .gitignore               # Excludes generated results and caches
└── tests/
    └── test_bet_rotor_mu_sweep.py # Comprehensive unit test suite
```

---

## Documentation

Full derivations, coordinate conventions, and analytical integrals are documented in detail in [zBET-documentation.md](zBET-documentation.md).
