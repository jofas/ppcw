FFLAGS=-g -O0 -check uninit,bounds -no-vec -fpp #-Dtest
FC=ifort

all: main
	make -C diff-output/
	make -C old/

clean:
	make -C diff-output/ clean
	make -C old/ clean

main: src/main.f90
	mkdir -p bin
	$(FC) $^ -o bin/$@ $(FFLAGS)

