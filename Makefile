FFLAGS=-g -O0 -check uninit,bounds -no-vec -fpp
FC=ifort

all: bench test test_all
	@mkdir -p bin
	@make -C diff-output/
	@make -C old/

clean:
	@make -C diff-output/ clean
	@make -C old/ clean

bench: src/main.f90
	$(FC) $^ -o bin/$@ $(FFLAGS)

test: src/main.f90
	$(FC) $^ -o bin/$@ $(FFLAGS) -Dtest

test_all: src/main.f90
	$(FC) $^ -o bin/$@ $(FFLAGS) -Dtest_all
