# AGENTS.md

## Repository Overview
`zBET` is a fast Python-based Blade Element Theory (BET) solver for helicopter and rotary-wing rotors in hover and forward flight across advance ratio sweeps ($\mu$). It is specifically designed for rapid conceptual sizing, trade studies, and parametric sweeps.

## Key Files
- `zBET.py`: Main solver containing rotor configuration, analytical BET equations, numerical integration, and export functions.
- `zBET-documentation.md`: Complete theoretical derivations, coordinate frame diagrams, sign conventions, and analytical formulas.
- `README.md`: High-level project summary, coordinate conventions, and quickstart guide.
- `tests/test_bet_rotor_mu_sweep.py`: Comprehensive standalone test suite (22 unit tests).
- `.gitignore`: Ensures all generated artifacts (`outputs/`, CSVs, PNGs, temporary files) remain untracked.

## Developer & Agent Guidelines
1. **Language**: All code comments, docstrings, terminal messages, and documentation must be written in **English**.
2. **Methodology**: Maintain the fast analytical BET architecture (closed-form weighted moments and global inflow gradient models) without introducing heavy iterative multi-annulus loops.
3. **Testing**: Always run `pytest` to ensure all 22 tests pass before proposing or committing changes.
4. **Execution**: Run `python zBET.py` to execute the forward flight sweep and generate outputs locally.
5. **Outputs**: Generated plots and CSV files belong in `outputs/` and must never be committed to version control.
6. **Simplicity**: Maintain a clean, minimalist repository structure. Avoid adding unnecessary root files or dependencies.
