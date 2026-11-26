# Contributing

## Setup

```bash
python -m venv .venv && source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -e ".[all,dev]"
pytest && ruff check .
```

## Rules

1. **Physics changes need a test against an independent reference.** Examples include brute
   force, a textbook case, or published data. "It runs" is not a test.
2. **Never weaken the calibration test** (`test_uncertainties_are_calibrated`) to make a
   change pass.
3. **Real-event results go in `docs/validation.md`**, including failures.
4. Keep the core dependency-free beyond NumPy and SciPy. Heavy dependencies go in extras.

## Most wanted

- Hand-reviewed sonic-boom picks for past re-entries (CSV format: see `src/falltrace/io.py`)
- Atmospheric acoustics (ray tracing, winds)
- Additional validation events with independent ground truth
