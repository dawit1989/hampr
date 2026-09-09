# Numerical Validation

Reproduces the **stage-by-stage** numerical comparison of HAMPR vs. the original
passive-radar benchmark on VEGA datasets 494–498.

## Files

| File | Purpose |
|------|---------|
| `bench_dump.cpp`   | Dumps every DSP stage from the **unmodified** benchmark (`.o` files linked). |
| `compare_dumps.cpp`| Compares two dump trees: exact-match count, max|Δ|, RMS, ULP. |
| `run_validation.sh`| End-to-end: build all tools, run both dumps, emit `comparison_report.md`. |
| `readres.cpp`      | Reads a `08_results.bin` dump (utility). |
| `Makefile`         | Builds `bench_dump` and `compare_dumps` manually. |

`tests/numerical_dump.cpp` is the HAMPR-side dump driver (built by the `numerical_dump`
CMake target). `tests/test_numerical_regression.cpp` is the final-output regression test.

## Requirements

- WSL2 / Linux with g++ 13, cmake, system FFTW3.
- Benchmark at `Code/passiveradar` (pre-built `*.o` + bundled `fftw/`).
- VEGA data at `Code/passiveradar_data/dataset` (datasets 494–498).

## Run

```
cd validation
./run_validation.sh        # full end-to-end run (~45s)
```

Or build tools manually:

```
make -C validation all
```

## Outputs (under `val_results/`, gitignored)

- `hampr/dump/<dataset>/<stage>.bin`  — HAMPR per-stage intermediates
- `bench/dump/<dataset>/<stage>.bin`  — benchmark per-stage intermediates
- `comparison_report.md`              — stage-by-stage metrics

See `docs/NUMERICAL_VALIDATION_REPORT.md` for the full analysis.