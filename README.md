# zBET

**zBET** is a fast, semi-empirical **Blade Element Theory (BET)** rotor solver for conceptual sizing, parametric sweeps, and preliminary rotor-performance studies. It combines closed-form radial moments for the main blade-element loads with global momentum theory for mean induced velocity and optional first-harmonic inflow models.

For the equations, assumptions, implementation mapping, and literature cross-check against Wayne Johnson and J. Gordon Leishman, see [zBET Documentation](zBET-documentation.md).

---

## Aerodynamic Models

zBET exposes the following aerodynamic selectors:

- `INDUCED_TORQUE_MODEL = "analytical_bet"`: direct BET torque integral for $C_{Qi}$.
- `INDUCED_TORQUE_MODEL = "energy_balance"`: infer $C_{Qi}$ from
  $K_{\mathrm{ind}}\lambda_iC_T+\mu_zC_T-\mu C_{Hi}$ without evaluating the direct induced-torque integral.
- `PROFILE_DRAG_MODEL = "analytical_tangential"`
- `PROFILE_DRAG_MODEL = "analytical_vectorial"`
- `PROFILE_DRAG_MODEL = "numerical_vectorial"`

Defaults:

```python
INDUCED_TORQUE_MODEL = "energy_balance"
PROFILE_DRAG_MODEL = "numerical_vectorial"
K_IND = 1.15
```

$K_{\mathrm{ind}}$ is used in energy-balance quantities: the energy-balance torque mode, `CPair`, and hover FoM. It is not inserted into the direct BET torque integral.

The principal outputs are:

| Quantity | Formulation |
| --- | --- |
| $C_T$ | BET |
| $C_{Hi}$ | BET |
| $C_{Qi}$ | selected direct-BET or energy-balance closure |
| $C_{H0}$, $C_{Q0}$ | selected profile-drag model |
| $C_Q$ | $C_{Qi}+C_{Q0}$ |
| $C_{P0,\mathrm{air}}$ | $C_{Q0}+\mu C_{H0}$ |
| $C_{Pair}$ | $K_{\mathrm{ind}}\lambda_iC_T+\mu_zC_T+C_{P0,\mathrm{air}}$ |
| FoM | hover energy balance with $K_{\mathrm{ind}}$ |


---

## BET vs. BEMT

zBET does not solve a separate local momentum balance at every annulus. Instead, it:

1. evaluates blade-element loads using analytical radial moments;
2. solves one global momentum-theory equation for mean induced velocity;
3. applies optional first-harmonic inflow gradients; and
4. evaluates profile drag with the selected tangential, vectorial analytical, or numerical-vectorial model.

This keeps execution fast while preserving the main physics needed for conceptual rotor studies.

---

## Coordinate System and Flight Conventions

- **Hub axes**
  - $+x$: forward.
  - $+y$: right / starboard.
  - $+z$: downward through the rotor disk.
- **Rotor rotation**
  - Viewed from above, the rotor turns counter-clockwise.
  - The advancing blade is on the right at $\psi=90^\circ$.
  - The local tangential velocity is $u_T=x+\mu\sin\psi$.
- **Thrust and drag**
  - Positive thrust acts upward, along $-z$.
  - Positive $C_H$ is rotor drag acting aft, along $-x$.
- **Shaft torque**
  - $C_Q>0$ is the positive magnitude of shaft torque required to power the rotor.
  - Shaft power is $P=Q\Omega$.
- **Axial flow**
  - $\lambda=\mu_z+\lambda_i$.
  - $\mu_z>0$: imposed relative flow is downward through the disk.
  - $\mu_z<0$: imposed relative flow is upward through the disk.
  - With `AXIAL_FLOW = "alpha"`, $\mu_z=-\mu\tan\alpha$.

---

## Features

- **Inflow models**
  - Uniform
  - Coleman simple
  - Coleman-Feingold / NDARC
  - Drees
- **Blade geometry**
  - Reference, physical geometric, and thrust-weighted solidity
  - Root cutout
  - Constant chord or linear taper
- **Pitch and hover trim**
  - Constant collective or linear twist
  - Collective trim to target $C_T$ or thrust
  - RPM trim to dimensional thrust
- **Engineering corrections**
  - Fixed or Sissingh-style effective tip-loss radius
  - Optional Prandtl-Glauert lift-slope correction
- **Outputs**
  - $C_T,C_Q,C_{Qi},C_{Q0},C_H,C_{Hi},C_{H0},C_Y,C_{Mx},C_{My}$
  - $C_{Pair}=K_{\mathrm{ind}}\lambda_iC_T+\mu_zC_T+C_{P0,\mathrm{air}}$
  - $C_{P0,\mathrm{air}}=C_{Q0}+\mu C_{H0}$
  - Effective rotor $L/D$
  - Hover figure of merit
  - Dimensional thrust and power

---

## Quick Start

### Requirements

- Python 3.9+
- `numpy`
- `pandas`
- `matplotlib`
- `pytest`

```bash
pip install numpy pandas matplotlib pytest
```

### Run the solver

```bash
python zBET.py
```

Generated CSV files and plots are written to `outputs/`.

### Run the test suite

```bash
pytest
```

---

## Project Structure

```text
zBET/
├── zBET.py
├── zBET-documentation.md
├── README.md
├── AGENTS.md
├── .gitignore
└── tests/
    └── test_bet_rotor_mu_sweep.py
```

---

## Documentation

The implementation-matched theory, sign conventions, mode definitions, assumptions, and literature references are in [zBET-documentation.md](zBET-documentation.md).
