# HAMPR Phase 8 — Multi-Static Processing

> Continuing from Phase 7 (Multi-Antenna / Multi-Receiver Scaling). Phase 6 (SDR) skipped.

## 1. Goal & Success Criteria

**Goal:** Add multi-static processing — FDoA estimation, hybrid TDoA/FDoA localization, GDOP quality metrics, iterative refinement.

**Success Criteria:**

| # | Criterion |
|---|-----------|
| 1 | FDoA estimation accuracy: within 10 Hz of true frequency offset |
| 2 | Hybrid TDoA/FDoA localization within 200 m of ground truth |
| 3 | GDOP computed for arbitrary receiver geometries |
| 4 | Iterative localization converges in <=3 iterations |
| 5 | All new code covered by tests (>=90% of new lines) |
| 6 | Existing tests continue to pass unchanged |

## 2. Current State (Phase 7 baseline)

- MultiStaticCorrelator — cross-correlation TDoA between site pairs
- TDoALocalizer — TDoA geolocation (single linearization step)
- MultiSitePipeline — orchestrates per-site processing + fusion

**Limitations:** Single linearization (no iteration), no FDoA, no GDOP, no hybrid fusion.

## 3. Architecture

Phase 8 extends Phase 7:
- MultiStaticCorrelator: add FDoA estimation (phase slope of cross-spectrum)
- New MultiStaticFusion class: hybrid TDoA/FDoA localization, GDOP, iterative refine

## 4. Implementation

### 4.1 FDoA Estimation
Split signal into two halves, compute cross-spectrum for each half.
FDoA = delta_phase / (2*pi * T_half), where T_half = N/(2*fs).

### 4.2 MultiStaticFusion class
- estimate_fdoa(ref_a, ref_b, fs) -> double (Hz)
- hybrid_localize(batch, track) -> GeoResult
- compute_gdop(site_positions, target_enu) -> double
- iterative_refine(initial, measurements, max_iter=3) -> refined position

### 4.3 Hybrid TDoA/FDoA Localization (Gauss-Newton)
1. Start with TDoA-only solution
2. Add FDoA equations
3. Linearize, solve via least-squares
4. Iterate until convergence

### 4.4 GDOP
GDOP = sqrt(trace((H^T H)^-1))

## 5. Test Plan

| Test | Description |
|------|-------------|
| Test 1 | FDoA estimation: known offset, error < 10 Hz |
| Test 2 | GDOP: known geometry, expected value |
| Test 3 | Hybrid TDoA/FDoA localization: error < 200 m |
| Test 4 | Iterative refinement: converges in <=3 iterations |
| Test 5 | Integration with MultiSitePipeline |

## 6. File Changes

New: include/hampr/dsp/multi_static_fusion.hpp
New: src/dsp/multi_static_fusion.cpp
New: tests/test_phase8.cpp
Modified: include/hampr/dsp/multi_static_correlator.hpp (FDoA methods)
Modified: src/dsp/multi_static_correlator.cpp (FDoA implementation)
Modified: CMakeLists.txt (new sources + test)
Modified: include/hampr/core/multi_site_types.hpp (FusionQuality struct)
Modified: include/hampr/dsp/multi_site_pipeline.hpp (MultiStaticFusion member)
Modified: src/dsp/multi_site_pipeline.cpp (integrate FDoA)

## 7. Implementation Order

| Step | Component | Duration |
|------|-----------|----------|
| 1 | FDoA estimation in MultiStaticCorrelator | 1 day |
| 2 | MultiStaticFusion class | 1 day |
| 3 | Integration with MultiSitePipeline | 0.5 day |
| 4 | Tests | 1 day |
| 5 | Documentation | 0.5 day |
