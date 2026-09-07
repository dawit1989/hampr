# HAMPR — Hardware-Accelerated Multi-Static Passive Radar

A modular, reproducible, hardware-accelerated passive-radar signal-processing platform based on the MilSpec Passive Radar benchmark and pyAPRiL implementation.

## Overview

HAMPR processes passive radar IQ data through a six-stage pipeline:

1. **Channel Isolation** — Separate reference and surveillance channels
2. **Beamforming** — Wiener/MVDR-style array processing
3. **Clutter Cancellation** — Time-domain Wiener-SMI-MRE filter
4. **Detection** — Overlap-and-Save cross-correlation with Hann windowing
5. **Target Finding** — Peak search in Range-Doppler map
6. **DOA Estimation** — MUSIC with eigen decomposition

### Phase 7: Multi-Antenna / Multi-Receiver Scaling

HAMPR supports multi-site, multi-antenna passive radar networks:

- **MultiSiteDataSource** — Loads IQ data from multiple receiver sites with JSON metadata
- **Synchronizer** — Clock offset correction and gain/phase calibration
- **ArrayGeometry** — Generalized scanning vectors for ULA, UCA, URA, and arbitrary arrays
- **MultiSitePipeline** — Per-site processing with cross-site fusion and TDoA localization
- **ResultFusion** — Bearing intersection and SINR-weighted combination
- **TDoALocalizer** — TDoA cross-correlation with WGS84/ECEF geolocation
- **MultiSiteGraph** — Streaming scheduler for multi-site batch processing

Single-site mode is fully backward-compatible: when `Config::multi_site_receivers` is empty, output is identical to the single-site `Pipeline`.

## Requirements

- C++17 compiler (GCC, Clang, MSVC)
- CMake 3.16+
- FFTW3 (`libfftw3-dev` on Debian/Ubuntu)
- Eigen3 (`libeigen3-dev` on Debian/Ubuntu)
- VEGA dataset for regression testing

### Optional (GPU acceleration)

- **SYCL backend**: A SYCL compiler (DPC++, hipSYCL, or ComputeCpp)
- **OpenCL backend**: OpenCL headers and runtime

## Build

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build

# With tests
cmake -B build -DHAMPR_BUILD_TESTS=ON
cmake --build build
cd build && ctest
```

### GPU Backends

```bash
# SYCL (CUDA/HIP/Metal via single-source C++)
cmake -B build -DHAMPR_ENABLE_GPU=SYCL
cmake --build build

# OpenCL fallback
cmake -B build -DHAMPR_ENABLE_GPU=OPENCL
cmake --build build

# CPU-only (default)
cmake -B build -DHAMPR_ENABLE_GPU=OFF
cmake --build build
```

## Usage

```bash
# Run the benchmark
./build/bin/benchmark_app

# Run tests
cd build && ctest
```

## Project Structure

```
hampr/
  CMakeLists.txt              # Build system
  PLAN.md                     # Engineering plan (all phases)
  PLAN_PHASE7.md              # Phase 7 detailed implementation plan
  include/hampr/
    core/
      types.hpp               # Type aliases (complex, array, mat, etc.)
      config.hpp              # Config struct (with multi-site fields)
      multi_site_types.hpp    # Multi-site data types
      flat_buffer.hpp         # Flat buffer data model (GPU-compatible)
    accel/
      fft_backend.hpp         # FFT backend interface
      linalg_backend.hpp      # LinAlg backend interface
      sycl_fft_backend.hpp    # SYCL FFT backend
      sycl_linalg_backend.hpp # SYCL LinAlg backend
      opencl_fft_backend.hpp  # OpenCL fallback FFT backend
    dsp/                      # DSP stage classes + multi-site components
      array_geometry.hpp      # ULA, UCA, URA, arbitrary array geometries
      synchronizer.hpp        # Clock sync and calibration
      fusion.hpp              # Cross-site bearing fusion
      localization.hpp        # TDoA geolocation
      multi_site_pipeline.hpp # Multi-site pipeline orchestrator
    io/                       # DataSource/DataSink abstractions
      multi_site_data_source.hpp # Multi-site data loading
    utils/                    # Math, FFTW, LinAlg, JSON utilities
    pipeline/                 # Pipeline schedulers
      multi_site_graph.hpp    # Multi-site streaming scheduler
  src/dsp/                    # DSP stage + multi-site implementations
  src/accel/                  # Backend implementations
  src/io/                     # I/O implementations
  src/pipeline/               # Pipeline scheduler implementations
  tests/                      # Unit and regression tests
  examples/                   # Benchmark application
```

## Backend Architecture

HAMPR uses backend abstraction for FFT and linear algebra operations:

```
DataSource -> [Pipeline: ChannelIsolator -> Beamformer -> ClutterCanceler
               -> Detector -> MetricExtractor -> DOAEstimator] -> DataSink
                   |
               FFTBackend / LinAlgBackend
```

**CPU backends**: FFTW (FFT), Eigen (LinAlg)
**GPU backends**: SYCL (single-source C++, multi-vendor), OpenCL (fallback)

Backend selection is compile-time via `HAMPR_GPU_BACKEND` preprocessor macro:
- `OFF` (default): CPU backends only
- `SYCL`: SYCL GPU backends
- `OPENCL`: OpenCL fallback backend

## Flat Buffer Data Model

`FlatBuffer` provides contiguous row-major storage for `complex` data, enabling:
- Zero-copy GPU memory transfer
- Improved CPU cache locality
- Conversion to/from `mat`/`array` via `from_mat()`/`to_mat()`

## Multi-Site Architecture

```
MultiSiteDataSource -> Synchronizer -> MultiSitePipeline -> MultiSiteGraph
       |                |               |
       |                |               +-> [Pipeline x N] -> ResultFusion -> TDoALocalizer
       |                |                   (per-site)       |
       |                |                                     -> MultiSiteResult
       |                |
       |                +-> clock_offset, gain, phase calibration
       |
       +-> JSON metadata (receiver positions, array geometry)
```

## Testing

| Test | Description |
|------|-------------|
| `test_core` | Math utility functions |
| `test_dsp` | DSP stage tests |
| `test_dsp_extended` | Extended DSP tests including FFT round-trip |
| `test_regression` | Full pipeline on VEGA dataset |
| `test_accel` | Backend consistency and flat buffer tests |
| `test_streaming` | Phase 5 streaming (BufferPool, StreamingDataSource) |
| `test_multi_site` | Phase 7 multi-site (ArrayGeometry, Synchronizer, Fusion, TDoA) |

## License

This project is based on the MilSpec Passive Radar benchmark and pyAPRiL implementation.
