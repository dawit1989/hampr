# HAMPR Phase 7 — Multi-Antenna / Multi-Receiver Scaling

> **Skipped Phase 6 (SDR Integration)** per user direction. This phase jumps
> directly to multi-receiver scaling, building on the Phase 5 streaming
> infrastructure and Phase 4 backend abstraction.

---

## 1. Goal & Success Criteria

**Goal:** Generalize HAMPR from a single 4-channel receiver site to a
multi-site, multi-antenna passive radar network. Add receiver metadata,
inter-site synchronization, generalized array geometries, cross-site DOA
fusion, and TDoA-based geolocation.

**Success Criteria:**

| # | Criterion |
|---|-----------|
| 1 | Process IQ data from N ≥ 2 receiver sites simultaneously |
| 2 | Single-site mode produces output identical to current `Pipeline` (within 1e-12) |
| 3 | Inter-site clock synchronization corrects offsets up to ±100 µs |
| 4 | Array geometry is configurable: ULA, UCA, URA, and arbitrary element positions |
| 5 | DOA bearings from multiple sites are fused via bearing intersection |
| 6 | TDoA localization produces geolocation within 500 m of ground truth |
| 7 | All new code covered by tests (≥ 90% of new lines) |
| 8 | CPU multi-site output matches reference within 1e-12 for synthetic data |

---

## 2. Current State Analysis

### Limitations Being Addressed

- **Fixed 4-channel assumption**: `ChannelIsolator` separates 1 ref + 3 surv
  channels with no metadata about which channels belong to which physical
  antenna or site.
- **No receiver metadata**: No information about site position, clock offset,
  antenna geometry, or per-receiver gain/phase calibration.
- **No synchronization**: All data is assumed to be perfectly time-aligned.
- **ULA-only scanning**: `DOAEstimator::gen_ula_scanning_vectors` hardcodes
  `phase = 2π · d · cos(θ)`. No support for circular or planar arrays.
- **Single-site pipeline**: `Pipeline::process` handles one `IQMatrix`.
  No fusion or localization across sites.
- **Nested-vector data model** (`mat`, `IQMatrix`) persists in DSP stages;
  `FlatBuffer` was introduced in Phase 4 but not yet adopted in `pipeline.cpp`,
  `detector.cpp`, or `clutter_canceler.cpp`.

### Architecture to Build On

- **Backend abstraction** (Phase 4): `FFTBackend`, `LinAlgBackend` with
  flat-buffer interfaces. Multi-site processing can use these for GPU offload.
- **Streaming infrastructure** (Phase 5): `BufferPool`, `StreamingDataSource`,
  `PipelineGraph`. Multi-site streaming reuses `BufferPool` per-site.
- **Type aliases** (Phase 4): `complex`, `array`, `mat`, `IQBuffer`, `IQMatrix`,
  `RDMap` centralized in `types.hpp`.
- **FlatBuffer** (Phase 4): contiguous 2D complex storage for GPU compatibility.

---

## 3. Data Model Changes

### 3.1 New Types (`include/hampr/core/multi_site_types.hpp`)

```cpp
struct AntennaElement {
    double x, y, z;       // Position in meters (local ENU frame)
};

struct ArrayGeometry {
    enum class Type { ULA, UCA, URA, ARBITRARY };
    Type type;
    std::vector<AntennaElement> elements;
    double element_spacing = 0.0;    // For ULA: inter-element spacing (λ units)
    double radius = 0.0;             // For UCA: array radius (λ units)
    int rows = 0, cols = 0;          // For URA: grid dimensions
};

struct ReceiverInfo {
    std::string id;
    double latitude = 0.0;       // WGS84 degrees
    double longitude = 0.0;
    double altitude = 0.0;       // meters MSL
    double clock_offset = 0.0;   // seconds, relative to master clock
    double gain = 0.0;           // dB, per-receiver calibration
    double phase_offset = 0.0;   // radians, per-receiver phase calibration
    ArrayGeometry array_geometry;
    int ref_channel_index = 0;
    int num_channels = 4;
};

struct MultiSiteData {
    std::vector<ReceiverInfo> receivers;
    std::vector<IQMatrix> iq_data;        // One IQMatrix per site
    double fs = 0.0;                      // Shared sampling rate
    double master_timestamp = 0.0;        // Reference time (seconds)
};

struct SynchronizedBatch {
    std::vector<ReceiverInfo> receivers;
    std::vector<IQMatrix> aligned_data;    // Time-aligned IQ data per site
    double fs = 0.0;
    int batch_samples = 0;
};

struct MultiSiteResult {
    std::vector<ProcessingResult> site_results;  // Per-site results
    double latitude = 0.0;
    double longitude = 0.0;
    double altitude = 0.0;
    double geolocation_error = 0.0;
    double fusion_time = 0.0;
    double total_time = 0.0;
};
```

### 3.2 Config Extension (`config.hpp`)

Additive fields, no breakage:

```cpp
struct Config {
    // ... existing fields ...
    std::vector<ReceiverInfo> multi_site_receivers;  // Empty = single-site mode
    bool enable_fusion = false;                       // Enable multi-site fusion
    bool enable_localization = false;                 // Enable TDoA/geolocation
    std::string localization_method = "tdoa";         // "tdoa" or "triangulation"
};
```

---

## 4. Component Design

### 4.1 MultiSiteDataSource (`include/hampr/io/multi_site_data_source.hpp`)

Extends `DataSource` to load data from multiple receiver sites.

```cpp
class MultiSiteDataSource : public DataSource {
public:
    // Load a single site's IQ data
    bool load_site(const std::string& site_id, const std::string& filename);

    // Load receiver metadata from a JSON or header file
    bool load_metadata(const std::string& metadata_file);

    // Get synchronized multi-site data
    MultiSiteData get_multi_site_data();

    // Streaming interface (per-site StreamingDataSource)
    bool open_site(const std::string& site_id, const std::string& filename);
    IQMatrix next_batch_site(const std::string& site_id, size_t batch_size);
    bool has_more_site(const std::string& site_id) const;
    void reset_all();
    void close_all();

    // Existing DataSource interface (single-site, for backwards compat)
    IQMatrix load(const std::string& filename, double& fs) override;
    TargetTrack load_track(const std::string& filename) override;

private:
    std::map<std::string, std::unique_ptr<StreamingDataSource>> site_sources_;
    std::vector<ReceiverInfo> receiver_infos_;
};
```

**Metadata format** (JSON):
```json
{
  "fs": 2500000.0,
  "receivers": [
    {
      "id": "site_A",
      "latitude": -37.8472,
      "longitude": 145.0451,
      "altitude": 12.0,
      "clock_offset": 0.0,
      "gain": -1.5,
      "phase_offset": 0.03,
      "array_geometry": {
        "type": "ULA",
        "num_elements": 4,
        "element_spacing": 0.528
      }
    },
    {
      "id": "site_B",
      "latitude": -37.8430,
      "longitude": 145.0480,
      "altitude": 15.0,
      "clock_offset": 0.000012,
      "gain": -2.0,
      "phase_offset": 0.01,
      "array_geometry": {
        "type": "ULA",
        "num_elements": 4,
        "element_spacing": 0.528
      }
    }
  ]
}
```

### 4.2 Synchronizer (`include/hampr/dsp/synchronizer.hpp`)

Corrects clock offsets and time-aligns data from multiple sites.

```cpp
class Synchronizer {
public:
    Synchronizer(double max_clock_error = 1e-4);

    // Align data from multiple sites to a common time reference
    SynchronizedBatch synchronize(const MultiSiteData& data);

    // Apply per-receiver gain and phase calibration
    void calibrate(MultiSiteData& data);

    // Interpolate to correct sub-sample clock offsets
    IQBuffer interpolate_shift(const IQBuffer& signal, double sample_offset);

private:
    double max_clock_error_;
    std::unique_ptr<FFTBackend> fft_backend_;
};
```

**Synchronization algorithm:**
1. Compute sample offset = `clock_offset * fs` for each site.
2. Integer part: shift samples by `floor(offset)`.
3. Fractional part: linear interpolation (or windowed sinc for better accuracy).
4. Apply per-receiver gain and phase: `sample *= gain * exp(j * phase_offset)`.
5. Return `SynchronizedBatch` with aligned data.

### 4.3 ArrayGeometry (`include/hampr/dsp/array_geometry.hpp`)

Generalized scanning vector generation for any array geometry.

```cpp
class ArrayGeometry {
public:
    ArrayGeometry() = default;
    explicit ArrayGeometry(const ReceiverInfo::ArrayGeometry& geom);

    // Generate scanning vectors for a set of look directions
    // For ULA: phase = 2π * d * cos(theta)
    // For UCA: phase = 2π * r * cos(theta - phi_i)
    // For URA: phase = 2π * (x*sin(theta)*cos(phi) + y*sin(theta)*sin(phi))
    mat gen_scanning_vectors(const std::vector<double>& thetas,
                             const std::vector<double>& phis = {}) const;

    // FlatBuffer variant for GPU path
    FlatBuffer gen_scanning_vectors_fb(const std::vector<double>& thetas,
                                       const std::vector<double>& phis = {}) const;

    int num_elements() const { return static_cast<int>(elements_.size()); }
    std::vector<AntennaElement> elements() const { return elements_; }

    // Compute array manifold for a given look direction
    array manaifold(double theta, double phi = 0.0) const;

private:
    std::vector<AntennaElement> elements_;
    ReceiverInfo::ArrayGeometry::Type type_;
    double element_spacing_;
    double radius_;
};
```

**Scanning vector formulas:**

| Geometry | Formula |
|----------|---------|
| ULA (1D) | `a_i = exp(j * 2π * d_i * cos(θ))` |
| UCA (2D) | `a_i = exp(j * 2π * r * cos(θ - φ_i))` |
| URA (2D) | `a_i = exp(j * 2π * (x_i * sin(θ) * cos(φ) + y_i * sin(θ) * sin(φ)))` |
| Arbitrary | `a_i = exp(j * 2π * (x_i * u_x + y_i * u_y + z_i * u_z))` where `u` is the look-direction unit vector |

### 4.4 MultiSitePipeline (`include/hampr/dsp/multi_site_pipeline.hpp`)

Orchestrates per-site processing and cross-site fusion.

```cpp
class MultiSitePipeline {
public:
    explicit MultiSitePipeline(const Config& config);

    // Process multi-site data: per-site processing → fusion → localization
    MultiSiteResult process(const MultiSiteData& data,
                            const TargetTrack& track);

    // Process a synchronized batch (streaming mode)
    MultiSiteResult process_batch(const SynchronizedBatch& batch,
                                  const TargetTrack& track);

    // Get per-site pipeline timing
    std::vector<double> site_processing_times() const;

private:
    Config config_;
    std::vector<std::unique_ptr<Pipeline>> site_pipelines_;
    std::unique_ptr<Synchronizer> synchronizer_;
    std::unique_ptr<ResultFusion> fusion_;
    std::unique_ptr<TDoALocalizer> localizer_;

    // Process a single site (delegates to Pipeline)
    ProcessingResult process_site(int site_idx,
                                  const IQMatrix& iq_data,
                                  const TargetTrackPoint& ref_track,
                                  double fs, int total_samples);
};
```

**Processing flow:**
1. For each site: run `ChannelIsolator`, `Beamformer`, `ClutterCanceler`,
   `Detector`, `MetricExtractor`, `DOAEstimator` (same as `Pipeline::process`)
2. Apply `ArrayGeometry` for DOA scanning (replacing hardcoded ULA)
3. Collect DOA bearings from all sites
4. If `enable_localization`:
   - TDoA: cross-correlate reference signals between site pairs, convert to
     hyperbolic lines, intersect for geolocation
   - Triangulation: intersect bearing lines from DOA estimates
5. If `enable_fusion`:
   - Weighted average of DOA estimates (weight = SINR)
   - Combined detection metric

### 4.5 ResultFusion (`include/hampr/dsp/fusion.hpp`)

Combines results from multiple receiver sites.

```cpp
class ResultFusion {
public:
    struct FusedTarget {
        double latitude;
        double longitude;
        double altitude;
        double confidence;
        std::vector<double> contributing_bearings;
        std::vector<std::string> contributing_sites;
    };

    FusedTarget fuse_bearings(const std::vector<ProcessingResult>& results,
                              const std::vector<ReceiverInfo>& receivers);

    double combined_sinr(const std::vector<ProcessingResult>& results);

private:
    // Iterative intersection of bearing lines (least-squares)
    std::pair<double, double> intersect_bearings(
        const std::vector<std::pair<double, double>>& site_positions,
        const std::vector<double>& bearings);
};
```

**Fusion algorithm:**
1. Collect (bearing, site_position) pairs from all sites.
2. For each pair of bearing lines, compute intersection point.
3. Use least-squares to find the point that minimizes sum of squared
   perpendicular distances to all bearing lines.
4. Confidence = function of bearing intersection angles (closer to 90° = better).

### 4.6 TDoALocalizer (`include/hampr/dsp/localization.hpp`)

TDoA-based geolocation using cross-correlation of reference signals.

```cpp
class TDoALocalizer {
public:
    struct GeoResult {
        double latitude;
        double longitude;
        double altitude;
        double residual_error;
    };

    GeoResult localize(const SynchronizedBatch& batch,
                       const TargetTrack& track);

private:
    // Cross-correlate reference channels between two sites to find time delay
    double estimate_tdoa(const IQBuffer& ref_a, const IQBuffer& ref_b, double fs);

    // Convert TDOA to hyperbola and intersect with site positions
    std::pair<double, double> tdoa_to_position(
        const std::vector<std::pair<double, double>>& site_positions,
        const std::vector<double>& tdoas, double c, double ref_lat, double ref_lon);

    // WGS84 → ECEF conversion
    void wgs84_to_ecef(double lat, double lon, double alt,
                        double& x, double& y, double& z);
    void ecef_to_wgs84(double x, double y, double z,
                        double& lat, double& lon, double& alt);
};
```

**TDoA algorithm:**
1. For each pair of sites (i, j): cross-correlate reference channels.
2. Find peak of cross-correlation → time delay τ_ij.
3. Convert τ_ij to distance difference: `Δd = c * τ_ij`.
4. Set up hyperbolic equations: `|P - S_i| - |P - S_j| = Δd`.
5. Solve via iterative least-squares (Newton-Raphson) or linear approximation.
6. Convert ECEF solution back to WGS84.

### 4.7 MultiSiteGraph (`include/hampr/pipeline/multi_site_graph.hpp`)

Streaming scheduler for multi-site processing, extending `PipelineGraph`.

```cpp
class MultiSiteGraph {
public:
    explicit MultiSiteGraph(const Config& config);

    std::vector<MultiSiteResult> process_stream(
        MultiSiteDataSource& source,
        const TargetTrack& track,
        size_t batch_size);

    // Async multi-site batch
    std::future<MultiSiteResult> process_batch_async(
        const SynchronizedBatch& batch,
        const TargetTrack& track);

    // Metrics
    double avg_latency_ms() const;
    double throughput_samples_per_sec() const;
    size_t batches_processed() const;

private:
    std::unique_ptr<MultiSitePipeline> pipeline_;
    std::unique_ptr<Synchronizer> synchronizer_;
    std::vector<std::future<MultiSiteResult>> futures_;
    std::vector<double> latencies_ms_;
    std::atomic<size_t> batches_processed_{0};
    std::atomic<int> total_samples_processed_{0};
    std::chrono::steady_clock::time_point start_time_;
    mutable std::mutex results_mtx_;
};
```

**Key difference from `PipelineGraph`:** Each batch pulls data from all
sites via `MultiSiteDataSource::next_batch_site`, then synchronizes, then
processes through `MultiSitePipeline`.

---

## 5. Implementation Steps

### Step 1: Multi-site data types (Day 1)
- Create `include/hampr/core/multi_site_types.hpp`
- Add `multi_site_receivers` and localization flags to `Config`
- Update `CMakeLists.txt` (header-only, no new sources)
- Tests: `test_multi_site.cpp` — Type construction, serialization

### Step 2: ArrayGeometry (Day 1–2)
- Create `include/hampr/dsp/array_geometry.hpp` and `src/dsp/array_geometry.cpp`
- Implement ULA, UCA, URA, and arbitrary scanning vector generation
- Unit tests: verify scanning vectors match known formulas
- Backwards compat: `Pipeline` uses `ArrayGeometry` internally for ULA mode

### Step 3: MultiSiteDataSource (Day 2)
- Create `include/hampr/io/multi_site_data_source.hpp` and `.cpp`
- Load multi-site IQ data from separate files + JSON metadata
- Streaming per-site interface using existing `StreamingDataSource`
- Tests: load two sites, verify alignment, error handling

### Step 4: Synchronizer (Day 2–3)
- Create `include/hampr/dsp/synchronizer.hpp` and `.cpp`
- Clock offset correction, sub-sample interpolation
- Gain/phase calibration
- Tests: synthetic clock offset correction within 1e-12

### Step 5: MultiSitePipeline core (Day 3–4)
- Create `include/hampr/dsp/multi_site_pipeline.hpp` and `.cpp`
- Per-site processing using existing `Pipeline` instances
- `ArrayGeometry` integration in DOA estimation
- Single-site mode: identical output to `Pipeline` within 1e-12
- Tests: single-site regression, multi-site processing

### Step 6: ResultFusion (Day 4)
- Create `include/hampr/dsp/fusion.hpp` and `.cpp`
- Bearing intersection, weighted SINR combination
- Tests: synthetic bearing intersection accuracy

### Step 7: TDoALocalizer (Day 4–5)
- Create `include/hampr/dsp/localization.hpp` and `.cpp`
- TDoA cross-correlation, hyperbola intersection
- WGS84↔ECEF conversion
- Tests: synthetic TDoA geolocation within 500 m

### Step 8: MultiSiteGraph streaming (Day 5)
- Create `include/hampr/pipeline/multi_site_graph.hpp` and `.cpp`
- Multi-site batch processing with async support
- Tests: multi-site streaming, metrics, async

### Step 9: Integration tests (Day 5)
- `test_multi_site.cpp` — full integration test with synthetic two-site data
- Backwards compatibility test: single-site mode matches `Pipeline`
- Numerical correctness: multi-site output verified against reference
- Run all existing tests to confirm no regression

### Step 10: Documentation (Day 5)
- Update `README.md` with multi-site usage
- Update `PLAN.md` with Phase 7 completion status

---

## 6. Test Plan

### New test file: `tests/test_multi_site.cpp`

| Test | Description |
|------|-------------|
| Test 1 | MultiSiteData construction and field access |
| Test 2 | ArrayGeometry ULA scanning vectors match reference formula |
| Test 3 | ArrayGeometry UCA scanning vectors match reference formula |
| Test 4 | ArrayGeometry URA scanning vectors match reference formula |
| Test 5 | ArrayGeometry arbitrary element positions |
| Test 6 | Synchronizer corrects integer sample offset |
| Test 7 | Synchronizer corrects fractional sample offset (`err < 1e-6`) |
| Test 8 | Synchronizer applies gain/phase calibration |
| Test 9 | MultiSitePipeline single-site mode matches `Pipeline` (`err < 1e-12`) |
| Test 10 | MultiSitePipeline two-site processing completes |
| Test 11 | ResultFusion bearing intersection accuracy (`error < 100 m`) |
| Test 12 | TDoALocalizer synthetic geolocation (`error < 500 m`) |
| Test 13 | MultiSiteGraph streaming batch processing |
| Test 14 | MultiSiteGraph async metrics interface |
| Test 15 | MultiSiteDataSource loads metadata + IQ data |

### Existing tests (no regression):
- `test_core`, `test_dsp`, `test_dsp_extended`, `test_regression`,
  `test_accel`, `test_streaming` — all must continue to pass unchanged.

### Numerical correctness strategy:
- Single-site mode: `MultiSitePipeline::process` with one site must produce
  `ProcessingResult` identical to `Pipeline::process` within 1e-12.
- TDoA: synthetic data with known site separation and target position;
  verify recovered position within 50% of baseline separation.
- Bearing fusion: synthetic DOA from 2+ sites; verify intersection point
  within 100 m of true position.

---

## 7. Assumptions

1. **Single frequency**: All receiver sites operate at the same carrier
   frequency and sampling rate.
2. **Post-cancellation phase**: Multi-site fusion happens after per-site
   clutter cancellation and detection, not during raw signal processing.
3. **File-based input**: No live SDR (Phase 6 skipped). Multi-site data is
   loaded from files with JSON metadata.
4. **WGS84 coordinates**: All site positions are in WGS84; conversions
   use standard ECEF transform.
5. **Speed of light**: TDoA uses `c = 299792458 m/s`.
6. **Phase centers**: Antenna element positions are in a local ENU frame
   relative to each receiver site's nominal position.
7. **Backwards compatibility**: When `Config::multi_site_receivers` is empty,
   `MultiSitePipeline` delegates to `Pipeline` and output is identical.
8. **No GPU dependency**: All multi-site code works on CPU. GPU backends
   (Phase 4) are used where `create_fft_backend()` / `create_linalg_backend()`
   are already called.
9. **Metadata format**: JSON is used for receiver metadata (simple, human-readable,
   no additional dependencies required — parsed manually or with a minimal parser).
10. **Clock model**: Linear clock offset (one scalar per site). No drift
    (rate) compensation — that is left for future Phase 7.5.
