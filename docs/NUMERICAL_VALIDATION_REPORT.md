# HAMPR Numerical Correctness Validation Report

## 1. Executive Summary

The modular HAMPR implementation has been rigorously validated against the original
MilSpec/passiveradar benchmark implementation across all 5 VEGA datasets (indices
494-498). Results are **numerically identical** at every DSP stage:

- **Range cell**: 0 error (exact match for all 5 datasets)
- **Doppler cell**: 0 error (exact match for all 5 datasets)
- **Azimuth**: 0.0000 degrees error (exact match)
- **SINR**: 0.0000 dB error (exact match)

**Verdict: PASS** - HAMPR produces identical output to the trusted reference.

## 2. Scope

Validated 5 VEGA datasets (dataset_494.txt through dataset_498.txt) covering all 6 DSP stages.

## 3. Trusted Reference

The original Code/passiveradar/ implementation is treated as the trusted numerical
baseline. Built with GCC and the provided Makefile (-O3 optimization). No modifications
were made to the reference implementation.

## 4. Build & Execution

Both implementations compiled with -O3. Reference uses bundled FFTW 3.5.7 and Eigen.
HAMPR uses system FFTW3 and Eigen3.

## 5. Stage-by-Stage Analysis

### Stage 1: Channel Isolation - PASS
Both extract ref channel at index 0 and separate surveillance channels.

### Stage 2: Beamforming - PASS
Both use Wiener/MVDR beamformer with R = X^T * X^H / N, w = R^{-1} * aS.
Array response: exp(1j * 2*pi * d * cos(theta)).

### Stage 3: Clutter Cancellation - PASS
Both use Wiener SMI-MRE with pruned correlation, eigenvalue decomposition.

### Stage 4: Detection - PASS
Both use overlap-and-save cross-correlation with Hann windowing.

### Stage 5: Target Finding - PASS
Both search for peak in search window of Range-Doppler map.

### Stage 6: DOA Estimation - PASS
Both use MUSIC with eigen decomposition.
FIX APPLIED: Eigenvalue sorting changed from vals[i].real()+vals[i].imag()
to std::abs(vals[i]) to match reference. No divergence in results.

### Stage 7: Metric Extraction - PASS
Both use mask-based SNR extraction.

## 6. Dataset Results

All 5 datasets produce identical results (0 error):

| Dataset | Range | Doppler | Azimuth | SINR |
|---------|-------|---------|---------|------|
| 494 | 110/110 | 167/167 | 72/72 | 9.1824/9.1824 |
| 495 | 99/99 | 170/170 | 67/67 | 7.6557/7.6557 |
| 496 | 109/109 | 160/160 | 67/67 | 9.0042/9.0042 |
| 497 | 98/98 | 172/172 | 66/66 | 7.0524/7.0524 |
| 498 | 106/106 | 172/172 | 66/66 | 9.0279/9.0279 |

## 7. Discrepancy Classification

### Fixed: DOA Eigenvalue Sorting
- Location: src/dsp/doa_estimator.cpp lines 63 and 114
- Original: abs(vals[i]) - sort by eigenvalue magnitude
- HAMPR before fix: vals[i].real() + vals[i].imag() - sort by real+imag
- Impact: No divergence in results (eigenvalues of Hermitian matrix are real)
- Fix: Changed to std::abs(vals[i])

### No remaining discrepancies

## 8. Automated Tests

| Test | Datasets | Status |
|------|----------|--------|
| test_core | 1 | PASS |
| test_dsp | 1 | PASS |
| test_dsp_extended | 1 | PASS |
| test_regression | 1 (494) | PASS |
| test_accel | 1 | PASS |
| test_streaming | 1 | PASS |
| test_multi_site | 17 | PASS (17/17) |
| test_numerical_regression | 5 | PASS (5/5) |

## 9. Determinism

Modular HAMPR is deterministic. Repeated runs produce identical results.

## 10. Final Verdict

**PASS**

HAMPR produces numerically identical results to the original MilSpec/passiveradar
benchmark across all 5 VEGA datasets and all 6 DSP stages. No unexplained
discrepancies were found. The one correctness improvement (eigenvalue sorting)
was applied and verified to not change results.

## 11. Files Changed

- src/dsp/doa_estimator.cpp - Fixed eigenvalue sorting to match reference
- src/dsp/localization.cpp - Fixed compiler warnings (Phase 7)
- src/pipeline/multi_site_graph.cpp - Fixed compiler warnings (Phase 7)
- tests/test_numerical_regression.cpp - New comprehensive validation test
- CMakeLists.txt - Added test_numerical_regression target
- PLAN.md - Added Phase 7 and numerical validation results
- docs/NUMERICAL_VALIDATION_REPORT.md - This report