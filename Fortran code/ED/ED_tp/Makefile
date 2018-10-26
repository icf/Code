FC = mpif90
PROG = lanh
FFLAGS=-DPOINTER_SIZE=8 -D MPI -fopenmp -I/home/xiaozhiyu/Documents/sprng2.0/include 
LIBS= /home/xiaozhiyu/Documents/sprng2.0/lib/libsprng.a \
      /home/xiaozhiyu/Documents/lapack-3.8.0/liblapack.a \
      /home/xiaozhiyu/Documents/lapack-3.8.0/librefblas.a\
      /usr/local/lib/libgmp.a 

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

.SUFFIXES: .o .for .F90
.F90.o : 
	$(FC) $(FFLAGS) -c  $<

.for.o : 
	$(FC) $(FFLAGS) -c  $<

.f.o : 
	$(FC) $(FFLAGS) -c  $<


.c.o : 
	$(CC) -O2  -c  $<

