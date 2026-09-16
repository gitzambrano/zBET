# zBET

**zBET** is a fast, semi-empirical rotor aerodynamic solver based on **Blade Element Theory (BET)** with analytical moment integration, coupled with global momentum theory for mean downwash. It is specifically designed for rapid conceptual sizing, parametric sweeps, and preliminary trade studies of helicopter and rotary-wing rotors.

> **BET vs. BEMT**: `zBET` evaluates blade section aerodynamics via closed-form radial-moment integrals and global inflow distributions rather than discretized, iterative multi-annulus BEMT strip-theory loops. This delivers execution times in milliseconds, ideal for conceptual design and optimizer sweeps.

For the complete theoretical derivations, coordinate frame diagrams, and analytical formulas, see the [zBET Documentation](zBET-documentation.md).

---

## Coordinate System & Flight Conventions

- **Hub Coordinate Axes**:
  - $x$-axis: Points **FORWARD** (aircraft nose / nominal flight direction).
  - $y$-axis: Points to the **RIGHT** (starboard side).
  - $z$-axis: Points **DOWNWARD** (through the bottom of the rotor disk).
- **Rotor Rotation**:
  - Viewed from **ABOVE** (looking down along $+z$): The blades rotate **COUNTER-CLOCKWISE (CCW)**.
  - Advancing blade is on the **STARBOARD / RIGHT** side ($\psi = 90^\circ$), where tangential speed is $u_T = x + \mu\sin\psi$.
  - Retreating blade is on the **PORT / LEFT** side ($\psi = 270^\circ$), where tangential speed is $u_T = x - \mu$.
- **Drive Torque & Fuselage Reaction**:
  - Aerodynamic blade drag resists rotation in the **CLOCKWISE (CW)** direction ($+z$ axis).
  - The drive shaft torque ($Q$) supplied by the engine acts in the **$-z$ direction** (CCW drive).
  - Rotor shaft torque coefficient $C_Q > 0$ is defined as the positive magnitude of torque required to power the rotor ($P = Q\Omega > 0$).
  - The reaction torque exerted by the rotor on the fuselage is **CLOCKWISE ($+z$ direction)** viewed from above, counteracted by the tail rotor.
- **Axial Inflow & Wind Direction ($\mu_z$)**:
  - Total axial inflow along the downward $+z$ axis is $\lambda = \mu_z + \lambda_i$, where induced downwash $\lambda_i \ge 0$ is always directed downward ($+z$).
  - **$\mu_z > 0$ (or vertical velocity $w > 0$)**: Relative wind flows **DOWNWARD** through the disk (oncoming wind coming from **ABOVE** the rotor, e.g. vertical climb).
  - **$\mu_z < 0$ (or vertical velocity $w < 0$)**: Relative wind flows **UPWARD** through the disk (oncoming wind coming from **BELOW** the rotor, e.g. vertical descent).
  - **Forward Flight with Forward Tilt ($\alpha > 0$)**: Relative wind comes from **BELOW** the disk ($\mu_z = -\mu \tan\alpha < 0$).

---

## Features

- **4 Inflow Models**:
  - Classic Uniform Inflow (Glauert, 1926)
  - Simple Coleman Inflow (Coleman et al., 1945)
  - NDARC Coleman-Feingold Inflow with lateral gradient $K_y$ (NASA/TP-2009-215402)
  - Drees Inflow with longitudinal and lateral gradients (Drees, 1949)
- **3 Solidity Options & Blade Taper**:
  - Reference solidity $\sigma_{\mathrm{ref}}$ (extrapolated to hub $r=0$)
  - True geometric physical solidity $\sigma_{\mathrm{geom}}$ (actual blade area / disk area)
  - Thrust-weighted equivalent solidity $\sigma_{\mathrm{thrust}}$ ($r^2$-weighted)
  - Linear blade chord taper ($c_{\mathrm{root}}$ to $c_{\mathrm{tip}}$)
- **Blade Pitch & Hover Trim Modes**:
  - Collective trim preserving total twist invariant ($\Delta\theta = \text{const}$)
  - RPM trim targeting rotor thrust in Newtons ($T$)
  - Constant collective pitch or linear twist ($\theta_{\mathrm{root}}$ to $\theta_{\mathrm{tip}}$)
- **Advanced Aerodynamic Corrections**:
  - Blade tip loss factor $B$ (fixed or Sissingh self-adjusting formula)
  - Prandtl-Glauert compressibility correction $a(M)$ on section lift curve slope
- **Comprehensive Performance Outputs**:
  - Non-dimensional coefficients: $C_T, C_Q, C_{Qi}, C_{Q0}, C_H, C_{Hi}, C_{H0}, C_Y, C_{Mx}, C_{My}$
  - Performance metrics: Total air power $C_{P,\mathrm{air}} = C_Q + \mu C_H$, effective rotor lift-to-drag $(L/D)_{\mathrm{eff}} = \mu C_T / C_{P,\mathrm{air}}$, and hover Figure of Merit ($FoM$)
  - Dimensional outputs: Thrust $T$ [N], Shaft Power $P$ [kW], and Total Air Power $P_{\mathrm{air}}$ [kW]

---

## Quick Start

### Requirements
- Python 3.9+
- `numpy`, `pandas`, `matplotlib`, `pytest`

```bash
pip install numpy pandas matplotlib pytest
```

### Running the Solver
Run the forward flight sweep with utility helicopter default parameters:
```bash
python zBET.py
```

Outputs are automatically generated in the `outputs/` directory:
- `outputs/zBET.csv` (consolidated dataset with all inflow models)
- `outputs/zBET_<model>.csv` (individual CSV per inflow model)
- 14 performance plots vs. $\mu$ for each inflow model

### Running Tests
Run the unit test suite:
```bash
pytest
```

---

## Project Structure
```text
zBET/
├── zBET.py                  # Core fast BET solver and advance ratio sweep engine
├── zBET-documentation.md    # Complete mathematical theory and analytical derivations
├── README.md                # Project overview and quickstart guide
├── AGENTS.md                # Minimalist guidelines for AI agents
├── .gitignore               # Excludes generated results, caches, and local files
└── tests/
    └── test_bet_rotor_mu_sweep.py # Comprehensive unit test suite
```

---

## Documentation

Full derivations, coordinate conventions, and analytical integrals are documented in [zBET-documentation.md](zBET-documentation.md).
