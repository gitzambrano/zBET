# zBET

**zBET** is a fast, semi-empirical **Blade Element Theory (BET)** rotor solver for conceptual sizing, parametric sweeps, and preliminary rotor-performance studies. It combines closed-form radial moments for the main blade-element loads with global momentum theory for mean induced velocity and optional first-harmonic inflow models.

For the equations, assumptions, implementation mapping, and literature cross-check against Wayne Johnson and J. Gordon Leishman, see [zBET Documentation](zBET-documentation.md).

---

## Aerodynamic Method

zBET uses fixed aerodynamic paths rather than selectable torque/power closures:

- $C_{Qi}$ is obtained from the BET torque expression.
- $C_{Q0}$ and $C_{H0}$ are obtained from direct vectorial radial × azimuthal integration of profile drag.
- $C_Q=C_{Qi}+C_{Q0}$.
- `CPair` is computed separately from an energy balance.

The tangential-only profile formulas are retained in the theory documentation only as a comparison showing what is lost when radial velocity is omitted. They are not used by the solver.

The principal loads are:

| Quantity | Formulation |
| --- | --- |
| $C_T$ | Analytical weighted-moment BET |
| $C_{Hi}$, $C_Y$ | Analytical weighted-moment BET |
| $C_{Mx}$, $C_{My}$ | Analytical weighted-moment BET |
| $C_{H0}$, $C_{Q0}$ | Direct vectorial profile-drag quadrature |
| $C_{Qi}$ | Direct BET torque expression |
| $C_Q$ | $C_{Qi}+C_{Q0}$ |
| $C_{Pair}$ | Energy balance |


---

## BET vs. BEMT

zBET does not solve a separate local momentum balance at every annulus. Instead, it:

1. evaluates blade-element loads using analytical radial moments;
2. solves one global momentum-theory equation for mean induced velocity;
3. applies optional first-harmonic inflow gradients; and
4. uses direct vectorial numerical quadrature for profile drag.

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
  - $C_T,C_{T0},C_Q,C_{Qi},C_{Q0},C_H,C_{Hi},C_{H0},C_Y,C_{Mx},C_{My}$
  - $C_{P,\mathrm{air}}=C_Q+\mu C_H-\mu_z C_T$
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
