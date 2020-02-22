FFLAGS=-g -O0 -check uninit,bounds -no-vec
FC=ifort

all: main

main: src/main.f90
	mkdir -p bin
	$(FC) $^ -o bin/$@ $(FFLAGS)
