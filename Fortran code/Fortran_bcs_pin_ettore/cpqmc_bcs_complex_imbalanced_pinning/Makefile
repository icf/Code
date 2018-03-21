FC = mpif90 -DMPI -O3
PROG = cpqmc
LIBS=-L/usr/local/lib -lgmp -L/home/evitali/sprng2.0/lib -lsprng -Wl,--no-as-needed -L${MKLPATH} -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lm
#FFLAGS=-fpp \
       -DPOINTER_SIZE=8 -I/sciclone/home00/mingpu/mylib/sprng2.0/include
FFLAGS= -cpp -DPOINTER_SIZE=8 -I/home/evitali/sprng2.0/include -DUSE_MKL -fopenmp -m64 -I${MKLINCLUDE}
#FFLAGS= -x f95-cpp-input -DPOINTER_SIZE=8 -I/sciclone/home00/ettore/sprng2.0/include
#FFLAGS= -Mpreprocess -DPOINTER_SIZE=8 -mcmodel=medium -I/sciclone/home00/ettore/sprng2.0/include 
OBJ=timing_module.o sprng_rndm.o io_module.o matrixcal.o caldet.o matrix_cpmc.o\
      params.o initial_end.o lattice_inform.o initial_kv.o initial_phiT.o matrix_dc.o\
      initial_phi.o modified_GS.o population.o adjustET.o update.o backup_phi.o\
      meas_sisj.o meas_ninj.o meas_pairing.o meas_cicj.o meas_htwo.o measure.o  data_manip.o file_name.o back_propagation.o\
      cpmc.o rcpmc.o step_qmc.o main.o measure_dynamics.o dynamics.o print_observables.o measure_energy.o\
      response_functions.o green_particles.o green_holes.o bcs_update.o bcs_green_pure.o


all: $(PROG)

rmdat:
	rm -f *.dat

clobber:
	rm -f *.o  *.mod *.out *.dat

clean:clobber
	rm -f cpqmc

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

