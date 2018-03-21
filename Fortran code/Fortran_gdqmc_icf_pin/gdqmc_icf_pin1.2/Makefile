FC = gfortran 

PROG = cpqmc-pinning

LIBS= /home/xiaozhiyu/Documents/sprng2.0/lib/libsprng.a \
      /home/xiaozhiyu/Documents/lapack-3.8.0/liblapack.a \
      /home/xiaozhiyu/Documents/lapack-3.8.0/librefblas.a\
      /usr/local/lib/libgmp.a 

FLAGS= -DPOINTER_SIZE=8

HEADS= -I/home/xiaozhiyu/Documents/sprng2.0/include

MODSF90= sprng_rndm.F90

MODSO=   sprng_rndm.o

# without MODSO
OBJ =   adjustET.o back_propagation.o backup_phi.o \
        caldet.o cpmc.o data_manip.o \
        file_name.o fortran_bug.o initial_end.o \
        initial_kv.o initial_phi.o initial_phiT.o \
        io_module.o lattice_inform.o  main.o \
        matrixcal.o matrix_cpmc.o matrix_dc.o \
        meas_cicj.o meas_htwo.o meas_sisj.o \
        measure.o modified_GS.o params.o \
        population.o rcpmc.o \
        step_qmc.o timing_module.o update.o

#out MODSF90
FUNS=   adjustET.F90 back_propagation.F90 backup_phi.F90 \
        caldet.F90 cpmc.F90 data_manip.F90 \
        file_name.F90 fortran_bug.f90 initial_end.F90 \
        initial_kv.F90 initial_phi.F90 initial_phiT.F90 \
        io_module.F90 lattice_inform.F90  main.F90 \
        matrixcal.F90 matrix_cpmc.F90 matrix_dc.F90 \
        meas_cicj.F90 meas_htwo.F90 meas_sisj.F90 \
        measure.F90 modified_GS.F90 params.F90 \
        population.F90 rcpmc.F90 \
        step_qmc.F90 timing_module.F90 update.F90

$(PROG): $(OBJ) 
	$(FC) $(FFLAGS) $(HEADS) -o $(PROG) $(OBJ) $(MODSO) $(LIBS)

$(OBJ): $(MODSO)
	-$(FC) -c $(FFLAGS) $(HEADS) $(FUNS) $(LIBS)
	$(FC) -c $(FFLAGS) $(HEADS) $(FUNS) $(LIBS)

$(MODSO): $(MODSF90)
	$(FC) -c $(FFLAGS) $(HEADS) $(MODSF90) $(LIBS)

	

.PHONY : dataclean refrash realclean

dateclean:
	rm -f *.dat

refrash:
	rm -f *.mod
	$(PROG)

realclean:
	rm -f $(PROG) *.mod *.dat *.out *.o







