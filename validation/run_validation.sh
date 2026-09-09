#!/usr/bin/env bash
# run_validation.sh - End-to-end numerical validation of HAMPR vs. the benchmark.
#
# Builds the HAMPR dump driver + benchmark dump driver + comparator, runs BOTH
# implementations on VEGA datasets 494-498, and writes val_results/comparison_report.md.
#
# Requirements: WSL2 with g++ 13, cmake, system FFTW3.
# Benchmark at ../passiveradar (pre-built *.o + bundled fftw/).
# VEGA data at ../passiveradar_data/dataset.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HAMPR_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PASSIVE_DIR="$(cd "$HAMPR_DIR/Code/passiveradar" && pwd)"
VAL_ROOT="$HAMPR_DIR/val_results"

echo "[validate] SCRIPT_DIR  = $SCRIPT_DIR"
echo "[validate] HAMPR_DIR   = $HAMPR_DIR"
echo "[validate] PASSIVE_DIR = $PASSIVE_DIR"

echo "[validate] (1/4) building HAMPR numerical_dump..."
cmake -S "$HAMPR_DIR" -B "$HAMPR_DIR/wsl_build" -Wno-dev >/tmp/hampr_cfg.log 2>&1
cmake --build "$HAMPR_DIR/wsl_build" --target numerical_dump >/tmp/hampr_build.log 2>&1
echo "[validate]   done."

echo "[validate] (2/4) building benchmark bench_dump..."
g++ -std=c++17 -O3 \
  -I "$PASSIVE_DIR/include" -I "$PASSIVE_DIR/eigen" -I "$PASSIVE_DIR/fftw" \
  "$SCRIPT_DIR/bench_dump.cpp" \
  "$PASSIVE_DIR/channelPreparation.o" "$PASSIVE_DIR/beamform.o" \
  "$PASSIVE_DIR/clutterCancellation.o" "$PASSIVE_DIR/detector.o" \
  "$PASSIVE_DIR/directionEstimation.o" "$PASSIVE_DIR/metricExtract.o" \
  "$PASSIVE_DIR/helper.o" \
  -L "$PASSIVE_DIR/fftw" -lfftw3 \
  -o "$SCRIPT_DIR/bench_dump"
echo "[validate]   done."

echo "[validate] (3/4) building compare_dumps..."
g++ -std=c++17 -O2 "$SCRIPT_DIR/compare_dumps.cpp" -o "$SCRIPT_DIR/compare_dumps"
echo "[validate]   done."

echo "[validate] (4/4) running dumps + comparison (this takes ~45s)..."
rm -rf "$VAL_ROOT"
"$HAMPR_DIR/wsl_build/numerical_dump" "$VAL_ROOT/hampr"
( cd "$PASSIVE_DIR" && LD_LIBRARY_PATH="$PASSIVE_DIR/fftw" "$SCRIPT_DIR/bench_dump" "$VAL_ROOT/bench" )
"$SCRIPT_DIR/compare_dumps" "$VAL_ROOT/hampr" "$VAL_ROOT/bench" > "$VAL_ROOT/comparison_report.md"

echo "[validate] report written to $VAL_ROOT/comparison_report.md"
tail -n 3 "$VAL_ROOT/comparison_report.md"