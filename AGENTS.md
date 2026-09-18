# AGENTS.md

## Repository Overview
`zBET` is a fast Python-based Blade Element Theory (BET) solver for helicopter and rotary-wing rotors in hover and forward flight across advance ratio sweeps ($\mu$). It is specifically designed for rapid conceptual sizing, trade studies, and parametric sweeps.

## Key Files
- `zBET.py`: Main solver containing rotor configuration, analytical BET equations, numerical integration, and export functions.
- `zBET-documentation.md`: Complete theoretical derivations, coordinate frame diagrams, sign conventions, and analytical formulas.
- `README.md`: High-level project summary, coordinate conventions, and quickstart guide.
- `tests/test_bet_rotor_mu_sweep.py`: Standalone regression and physics-consistency test suite.
- `.gitignore`: Ensures all generated artifacts (`outputs/`, CSVs, PNGs, temporary files) remain untracked.

## Developer & Agent Guidelines
1. **Language**: All code comments, docstrings, terminal messages, and documentation must be written in **English**.
2. **Methodology**: Maintain the fast analytical BET architecture (closed-form weighted moments and global inflow gradient models) without introducing heavy iterative multi-annulus loops. Preserve the explicit selectors `INDUCED_TORQUE_MODEL` (`analytical_bet` or `energy_balance`) and `PROFILE_DRAG_MODEL` (`analytical_tangential`, `analytical_vectorial`, or `numerical_vectorial`). In direct BET torque mode, do not insert `K_IND` into the torque integral. In energy-balance torque mode, use `K_IND`. Always compute `CPair` from the energy balance with `K_IND`: `CPair = K_IND*lambda_i*CT + mu_z*CT + CP0_air`, with `CP0_air = CQ0 + mu*CH0`. Do not add a separate `mu*CHi` term to `CPair`; in `energy_balance` mode it is already accounted for through `CQi = K_IND*lambda_i*CT + mu_z*CT - mu*CHi`. Compute hover FoM with `K_IND`.
3. **Testing**: Always run the full `pytest` suite before proposing or merging behavioral changes. Do not hard-code the collected test count in documentation.
4. **Execution**: Run `python zBET.py` to execute the forward flight sweep and generate outputs locally.
5. **Outputs**: Generated plots and CSV files belong in `outputs/` and must never be committed to version control.
6. **Simplicity**: Maintain a clean, minimalist repository structure. Avoid adding unnecessary root files or dependencies.
