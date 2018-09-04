FC = ifort 
PROG = lanh
FFLAGS=-fpp -DPOINTER_SIZE=8 -I/home/boruoshihao/sprng/sprng2.0/include\
       -i8 -I/opt/intel/mkl/include -openmp
LIBS= -L/home/boruoshihao/sprng/sprng2.0/lib -lsprng -lgmp\
      -L/opt/intel/lib/intel64 -limf\
      -L/opt/intel/mkl/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread -lm
OBJ = io_module.o timing_module.o sprng_rndm.o params.o initial_end.o set_binomial.o\
      lattice_inform.o filename.o set_matrix_updn.o lanczos.o  measure.o  measure_corr.o measure_Square.o measure_k.o\
      measure_cicj.o project.o main.o


all: $(PROG)

rmdat:
	rm -f *.dat

clobber:
	rm -f *.o  *.mod *.out *.dat

clean:clobber
	rm -f lanh

$(PROG): $(OBJ)
	$(FC) $(FFLAGS) -o $(PROG) $(OBJ) $(LIBS)

.SUFFIXES: .o .for .f90
.f90.o : 
	$(FC) $(FFLAGS) -c  $<

.for.o : 
	$(FC) $(FFLAGS) -c  $<

.f.o : 
	$(FC) $(FFLAGS) -c  $<


.c.o : 
	$(CC) -O2  -c  $<

