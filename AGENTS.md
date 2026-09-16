# AGENTS.md

## Repository Overview
`zBET` is a Python-based Blade Element Momentum Theory (BEMT) solver for helicopter rotors in hover and forward flight across advance ratio sweeps ($\mu$).

## Key Files
- `zBET.py`: Main solver containing rotor configuration, BEMT equations, numerical integration, and export functions.
- `zBET-documentation.md`: Complete theoretical derivations, coordinate conventions, and analytical formulas.
- `README.md`: High-level project summary and quickstart guide.
- `tests/test_bet_rotor_mu_sweep.py`: Comprehensive test suite (22 unit tests).
- `.gitignore`: Ensures all generated artifacts (`outputs/`, CSVs, PNGs, temporary files) remain untracked.

## Developer & Agent Guidelines
1. **Language**: All code comments, docstrings, terminal messages, and documentation must be written in **English**.
2. **Testing**: Always run `pytest` to ensure all tests pass before proposing or committing changes.
3. **Execution**: Run `python zBET.py` to execute the forward flight sweep and generate outputs locally.
4. **Outputs**: Generated plots and CSV files belong in `outputs/` and must never be committed to version control.
5. **Simplicity**: Maintain a clean, minimalist repository structure. Avoid adding unnecessary root files or dependencies.
