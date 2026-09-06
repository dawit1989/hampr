# HAMPR — Engineering Plan

## Results: First Implementation Slice (Completed)

### Changes Applied

1. **fftw_wrapper.hpp** — Fixed ifft() to normalize by N. Added cached FFT plan mechanism.
2. **clutter_canceler.cpp** — Restored conjugate() before fftshift() in pruned_corr().
3. **doa_estimator.cpp** — Fixed eigenvalue sorting to match original.
4. **CMakeLists.txt** — Split hampr_tests into test_core, test_dsp, test_regression.
5. **test_regression.cpp** — New regression test for VEGA index 494.

### Correctness

The fixed modular implementation produces identical output to the original passiveradar/ implementation across all 5 VEGA time indices.

### Performance

| Stage | Before | After | Speedup |
|-------|--------|-------|---------|
| Clutter Cancellation | 1.067s | 0.233s | 4.6x |
| DOA / MUSIC | 2.301s | 0.685s | 3.4x |
| Total | 3.538s | 1.007s | 3.5x |

### Acceptance Criteria

- Project builds with cmake
- All tests pass (ctest)
- FFT plan caching implemented
- CPU output matches reference

## Phase 4 Implementation Results

### Completed Phase 4 Work

1. **Type alias centralization** — Moved array and mat to types.hpp as canonical definitions.
2. **Flat buffer data model** — FlatBuffer class with row-major contiguous storage.
3. **Backend abstraction** — FFTBackend and LinAlgBackend interfaces with flat buffer methods.
4. **SYCL GPU backend** — Single-source C++ using DPC++/hipSYCL/ComputeCpp.
5. **OpenCL fallback backend** — OpenCL C++ binding FFT backend.
6. **GPU backend factory** — Compile-time backend selection via HAMPR_GPU_BACKEND.
7. **CMakeLists.txt** — HAMPR_ENABLE_GPU cache variable, test_accel target.
8. **Phase 4 tests** — 14 tests covering all backends and flat buffer methods.

### Remaining Phase 4 Work

- **Flat buffer migration of DSP stages** — Replace mat/array with FlatBuffer in doa_estimator.cpp, detector.cpp, clutter_canceler.cpp.

## A. Current Architecture

### Repository Layout

Code/
  hampr/ — Modular C++17 implementation
  passiveradar/ — Original monolithic implementation
  passiveradar_data/ — VEGA dataset (494-498)

### Pipeline Architecture

| Stage | Original | Modular |
|-------|----------|---------|
| 1. Channel Isolation | isolate_channels() | ChannelIsolator::isolate() |
| 2. Beamforming | beamform_surveillance() | Beamformer::beamform() |
| 3. Clutter Cancellation | Wiener_SMI_MRE() | ClutterCanceler::wiener_smi_mre() |
| 4. Detection | cc_detector_ons() | Detector::cc_detector_ons() |
| 5. Target Finding | find_target() | MetricExtractor::extract() |
| 6. DOA Estimation | MUSIC + Eigen | DOAEstimator::estimate() |

### Data Types

types.hpp: complex, IQBuffer, IQMatrix, RDMap, array, mat, TargetTrackPoint, TargetTrack, ProcessingResult.

### Build System

CMake-based build with FFTW3 and Eigen3. Options: HAMPR_ENABLE_GPU, HAMPR_BUILD_TESTS, HAMPR_BUILD_EXAMPLES.

## B. Architecture Assessment

### Strengths

- Clean separation of stages
- Abstract I/O interfaces
- Per-stage timing
- Preserved algorithmic structure
- Original implementation as reference

### Weaknesses and Technical Debt

- Data model uses vector<vector> (not GPU-compatible)
- Hardcoded angular resolution (180-degree scan)
- No streaming DataSource
- DSP stages call utilities directly, not through backends

### Scalability Limitations

- Fixed 4-channel assumption
- No receiver metadata or synchronization
- All data loaded into memory

## C. Performance Assessment

Primary bottlenecks: DOA/MUSIC (eigen decomposition + scanning), Clutter Cancellation (FFT), Cross-Correlation Detection (FFT). Memory movement from nested vectors.

## D. Correctness Assessment

Validated against original passiveradar/ across all 5 VEGA indices. Results match exactly. Existing tests: test_core, test_dsp, test_dsp_extended, test_regression.

## E. Target Architecture

Backend abstraction: FFTBackend (CPU: FFTW, GPU: SYCL/OpenCL), LinAlgBackend (CPU: Eigen, GPU: SYCL). FlatBuffer data model for GPU compatibility.

## F. Phased Roadmap

### Phase 1 — Establish the Reference Baseline

Ensure builds, tests pass, VEGA dataset runs, outputs match. Document baseline.

### Phase 2 — Strengthen Modular Architecture

Improve interfaces, configuration, data structures, error handling, testing, I/O separation.

### Phase 3 — Profiling and CPU Optimization

Profile pipeline. Focus on DOA/MUSIC, clutter cancellation, cross-correlation, FFT, memory movement.

### Phase 4 — Hardware Acceleration (Multi-Architecture)

SYCL primary + OpenCL fallback. Flat buffer data model. Backend abstraction. CPU output matches reference within 1e-12.

### Phase 5 — Streaming

Streaming DataSource, buffering, asynchronous execution, latency/throughput measurements.

### Phase 6 — SDR Integration

Live SDR input through DataSource interface. No SDR coupling in DSP.

### Phase 7 — Multi-Antenna / Multi-Receiver Scaling

Generalize beyond 4 channels. Receiver metadata and synchronization.

### Phase 8 — Multi-Static Processing

Multi-receiver measurement handling, fusion and localization.

## G. Priority Ranking

| Work Item | Priority | Rationale |
|-----------|----------|-----------|
| Centralize type aliases in types.hpp | P0 | Build-breaking compilation error |
| Backend abstraction (FFT/LinAlg) | P1 | Enables GPU acceleration |
| Flat buffer data model | P1 | Required for GPU |
| CPU backends (FFTW/Eigen) | P1 | GPU backend fallback |
| Regression test | P0 | Numerical validation gap |
| SYCL GPU backend | P1 | Multi-architecture acceleration |
| OpenCL fallback backend | P2 | Fallback when SYCL unavailable |
| Streaming DataSource | P1 | Required for SDR |
| SDR integration | P3 | Requires hardware |
| Multi-static processing | P3 | Future expansion |

## H. First Implementation Slice

Objective: Establish verified reference baseline. Files: fftw_wrapper.hpp, clutter_canceler.cpp, doa_estimator.cpp, test_regression.cpp, CMakeLists.txt.
