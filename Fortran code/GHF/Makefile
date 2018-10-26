FC = pgf90  
#FC = mpif90 -DMPI
PROG = sxhf
#LIBS=-L/share/apps/intel/lib/intel64 -limf\
     -L/share/apps/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lguide -lpthread -lm\
     -L/home/mingpu/sprng2.0/lib -lsprng -lgmp

LIBS=-L${ACML}/lib -lacml -L/storm/home00/mingpu/software/sprng2.0/lib -lsprng -lgmp  -lm -L/usr/local/pgi-14.3/linux86-64/14.3/libso -lpgmp
#FFLAGS=-fpp \
       -DPOINTER_SIZE=8 -I/home/mingpu/sprng2.0/include
FFLAGS= -Mpreprocess -DPOINTER_SIZE=8 -I/storm/home00/mingpu/software/sprng2.0/include
#FFLAGS=-shared-intel -mcmodel=medium -fpp \
       -DPOINTER_SIZE=8 -I/home/boruoshihao/sprng/sprng2.0/include
OBJ = timing_module.o io_module.o sprng_rndm.o matrixcal.o caldet.o matrix_cpmc.o\
      params.o initial_end.o lattice_inform.o x-hf.o main.o
#      caldet.o io_module.o\
      matrixcal.o matrix_cpmc.o lattice_inform.o initial_end.o\
      diagonalize_H.o sort2.o degen_coup.o initial_phiT.o initial_cpmc.o\
      update.o population.o adjustET.o modified_GS.o gs_subs.o\
      forw_meas.o back_propagation.o measure_coor.o file_name.o HF.o data_manip.o\
      measure_cicj.o main.o
#      free_meas.o free_bk_meas.o free_bb_meas.o\


all: $(PROG)

rmdat:
	rm -f *.dat

clobber:
	rm -f *.o  *.mod *.out *.dat

clean:clobber
	rm -f sxhf

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

