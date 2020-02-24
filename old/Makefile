SRC=MD.f util.f
OBJ=$(SRC:.f=.o)

FFLAGS= -g -O0  -check uninit,bounds -no-vec -fpp
FC=ifort

all: old_bench old_test

old_bench: $(OBJ)
	$(FC) $(FFLAGS)  -o ../bin/$@ control.f $(OBJ)

old_test: $(OBJ)
	$(FC) $(FFLAGS) -Dtest -o ../bin/$@ control.f $(OBJ)

clean:
	rm -f $(OBJ)
	rm -f control.o

$(OBJ): coord.inc Makefile

MD.o: MD.f
	$(FC) $(FFLAGS) -c MD.f

control.o: control.f
	$(FC) $(FFLAGS) -c control.f

util.o: util.f
	$(FC) $(FFLAGS) -c util.f




