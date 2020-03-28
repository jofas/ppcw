PROG=src/history.f90
#PROG=src/main.f90
FFLAGS_PROFILE=-fno-inline-functions
FFLAGS_BENCH=-qopt-report=5
#FFLAGS=-Ofast -fpp -qopenmp -vec-threshold0 -xCORE-AVX2 \
	-align array64byte -ipo -qopt-prefetch=5
FFLAGS=-O3 -fpp -no-vec
FC=ifort

all: test test_all bench profile
	@mkdir -p bin
	@make -C diff-output/
	@make -C old/

clean:
	rm bin/*
	@make -C diff-output/ clean
	@make -C old/ clean

bench: $(PROG)
	$(FC) $^ -o bin/$@ $(FFLAGS) $(FFLAGS_BENCH)

profile: $(PROG)
	$(FC) $^ -o bin/$@ $(FFLAGS) $(FFLAGS_PROFILE)

test: $(PROG)
	$(FC) $^ -o bin/$@ $(FFLAGS) -Dtest

test_all: $(PROG)
	$(FC) $^ -o bin/$@ $(FFLAGS) -Dtest_all
