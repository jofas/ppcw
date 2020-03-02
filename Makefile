FFLAGS=-O3 -no-vec -fpp -p
FC=ifort

all: test test_all bench
	@mkdir -p bin
	@make -C diff-output/
	@make -C old/

clean:
	@make -C diff-output/ clean
	@make -C old/ clean

bench: src/main.f90
	$(FC) $^ -o bin/$@ $(FFLAGS) -qopt-report=5

test: src/main.f90
	$(FC) $^ -o bin/$@ $(FFLAGS) -Dtest

test_all: src/main.f90
	$(FC) $^ -o bin/$@ $(FFLAGS) -Dtest_all
