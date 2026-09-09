# HAMPR Phase 6 - SDR Integration

> Re-visiting after Phase 8. Originally skipped to prioritize Phase 7/8;
> core DSP and streaming infrastructure is now complete.

## 1. Goal & Success Criteria

**Goal:** Add live SDR input through the DataSource interface. No SDR coupling
in the DSP layer.

**Success Criteria:**

| # | Criterion |
||---|-----------|
|| 1 | RadioSource abstract interface for SDR hardware backends |
|| 2 | FileRadioSource mock for testing (file-based replay) |
|| 3 | LiveDataSource implements DataSource, uses RadioSource |
|| 4 | MultiSiteDataSource can consume live SDR sources |
|| 5 | Tests pass with file-based SDR mock (no hardware required) |
|| 6 | DSP code unchanged - no SDR headers in dsp/ includes |

## 2. Current State

- DataSource (abstract): load(), load_track()
- StreamingDataSource: file-based streaming via BufferPool
- MultiSiteDataSource: multi-site file loading with JSON metadata
- Pipeline, MultiSitePipeline, MultiStaticFusion all use IQ data

**Gap:** No live SDR input. All data comes from files.

## 3. Architecture

FileRadioSource (mock) -> RadioSource (abstract) -> LiveDataSource (DataSource impl) -> DSP

No SDR coupling in DSP layer.

## 4. Implementation

### 4.1 RadioSource (abstract)
Interface for SDR hardware:
- open(device), start(), stop(), read(), set_frequency(), set_gain(),
  set_sample_rate(), sampling_rate(), num_channels(), is_open()

### 4.2 FileRadioSource (mock)
File-based replay of IQ data for testing.

### 4.3 LiveDataSource
Implements DataSource, uses RadioSource via BufferPool.

## 5. Test Plan

| Test | Description |
||------|-------------|
|| Test 1 | FileRadioSource loads IQ data from file |
|| Test 2 | LiveDataSource streams data via BufferPool |
|| Test 3 | RadioSource interface (open/start/stop/read) |
|| Test 4 | No SDR headers in dsp/ includes |
|| Test 5 | MultiSiteDataSource with live sources |

## 6. File Changes

New: include/hampr/io/radio_source.hpp
New: src/io/radio_source.cpp
New: include/hampr/io/live_data_source.hpp
New: src/io/live_data_source.cpp
New: tests/test_phase6.cpp
Modified: CMakeLists.txt
Modified: README.md

## 7. Implementation Order

| Step | Component | Duration |
||------|-----------|----------|
|| 1 | RadioSource + FileRadioSource | 1 day |
|| 2 | LiveDataSource | 0.5 day |
|| 3 | MultiSiteDataSource integration | 0.5 day |
|| 4 | Tests | 1 day |
|| 5 | Documentation | 0.5 day |
