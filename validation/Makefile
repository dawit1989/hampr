# Makefile for validation tools.
# Usage: make -C validation all
#
# Searches for the benchmark tree in Code/passiveradar first, then falls back
# to the original source tree in Ref/MilSpec Bench Passive Radar/passiveradar.
THIS_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
PASSIVE_DIR ?= $(shell for d in \
    $(THIS_DIR)/../Code/passiveradar \
    $(THIS_DIR)/../Ref/MilSpec Bench Passive Radar/passiveradar \
    /mnt/d/HAMPR/Code/passiveradar \
    /mnt/d/HAMPR/Ref/MilSpec Bench Passive Radar/passiveradar; do \
    if [ -d "$$d/include" ] && [ -f "$$d/beamform.o" ]; then echo "$$d"; break; fi; \
done)

.PHONY: all bench_dump compare_dumps readres clean check

check:
	@if [ -z "$(PASSIVE_DIR)" ]; then \
    echo "ERROR: benchmark tree not found (need include/ + *.o + fftw/ + eigen/)"; \
    echo "Looked in: Code/passiveradar, Ref/MilSpec Bench Passive Radar/passiveradar"; \
    exit 1; \
	fi

all: bench_dump compare_dumps readres

bench_dump: bench_dump.cpp
	$(MAKE) check
	g++ -std=c++17 -O3 -Wall -pedantic -I$(PASSIVE_DIR)/include -I$(PASSIVE_DIR)/eigen -I$(PASSIVE_DIR)/fftw \
		bench_dump.cpp $(PASSIVE_DIR)/channelPreparation.o $(PASSIVE_DIR)/beamform.o \
		$(PASSIVE_DIR)/clutterCancellation.o $(PASSIVE_DIR)/detector.o \
		$(PASSIVE_DIR)/directionEstimation.o $(PASSIVE_DIR)/metricExtract.o \
		$(PASSIVE_DIR)/helper.o \
		-L$(PASSIVE_DIR)/fftw -lfftw3 -o bench_dump

compare_dumps: compare_dumps.cpp
	g++ -std=c++17 -O2 compare_dumps.cpp -o compare_dumps

readres: readres.cpp
	g++ -std=c++17 -O2 readres.cpp -o readres

clean:
	rm -f bench_dump compare_dumps readres
