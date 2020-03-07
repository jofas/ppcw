FFLAGS_BENCH=-qopt-report=5 #-p
FFLAGS=-O3 -fpp -qopenmp -vec-threshold0 -xCORE-AVX2 \
	-align array64byte #-mkl
FC=ifort

all: test test_all bench
	@mkdir -p bin
	@make -C diff-output/
	@make -C old/

clean:
	@make -C diff-output/ clean
	@make -C old/ clean

bench: src/main.f90
	$(FC) $^ -o bin/$@ $(FFLAGS) $(FFLAGS_BENCH)

test: src/main.f90
	$(FC) $^ -o bin/$@ $(FFLAGS) -Dtest

test_all: src/main.f90
	$(FC) $^ -o bin/$@ $(FFLAGS) -Dtest_all
