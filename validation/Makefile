# Makefile for validation tools.
# Usage: make -C validation all
PASSIVE_DIR ?= $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/../Code/passiveradar)

.PHONY: all bench_dump compare_dumps readres clean

all: bench_dump compare_dumps

bench_dump: bench_dump.cpp
	g++ -std=c++17 -O3 -I$(PASSIVE_DIR)/include -I$(PASSIVE_DIR)/eigen -I$(PASSIVE_DIR)/fftw \
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
