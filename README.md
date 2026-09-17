# zBET

**zBET** is a fast, semi-empirical **Blade Element Theory (BET)** rotor solver for conceptual sizing, parametric sweeps, and preliminary rotor-performance studies. It combines closed-form radial moments for the main blade-element loads with global momentum theory for mean induced velocity and optional first-harmonic inflow models.

For the equations, assumptions, implementation mapping, and literature cross-check against Wayne Johnson and J. Gordon Leishman, see [zBET Documentation](zBET-documentation.md).

---

## Aerodynamic Model Selectors

zBET uses **two independent and explicit selectors**:

- `PROFILE_DRAG_MODEL`
  - `"analytical_bet"`: closed-form BET profile force and torque.
  - `"numerical_profile"`: radial × azimuthal Gauss-Legendre profile-drag quadrature.
- `INDUCED_TORQUE_MODEL`
  - `"analytical_bet"`: direct analytical BET induced torque.
  - `"energy_balance"`: induced shaft torque from the energy-balance closure with `K_IND`.

The default configuration is:

```python
PROFILE_DRAG_MODEL = "numerical_profile"
INDUCED_TORQUE_MODEL = "energy_balance"
```

The principal lift-induced loads remain analytical in both selector families:

| Quantity | Formulation |
| --- | --- |
| $C_T$ | Analytical weighted-moment BET |
| $C_{Hi}$, $C_Y$ | Analytical weighted-moment BET |
| $C_{Mx}$, $C_{My}$ | Analytical weighted-moment BET |
| $C_{H0}$, $C_{Q0}$ | Selected by `PROFILE_DRAG_MODEL` |
| $C_{Qi}$ | Selected by `INDUCED_TORQUE_MODEL` |
| $C_Q$ | $C_{Qi}+C_{Q0}$ |

The selectors are deliberately named after the physics they implement. zBET does not expose a generic “simple/complete” switch.

A fully integrated force-balance solver would evaluate the same local aerodynamic state consistently for all force and moment channels over radius and azimuth. zBET is intentionally lighter than that: its design target is fast, transparent conceptual analysis.
---

## BET vs. BEMT

zBET does not solve a separate local momentum balance at every annulus. Instead, it:

1. evaluates blade-element loads using analytical radial moments;
2. solves one global momentum-theory equation for mean induced velocity;
3. applies optional first-harmonic inflow gradients; and
4. uses numerical quadrature only where explicitly selected for profile drag.

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
  - $C_{P,\mathrm{air}}=C_Q+\mu C_H$
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
