FC = mpif90

PROG = cpqmc-pinning

LIBS= /sciclone/pscr/zxiao01/lib/sprng2.0/lib/libsprng.a \
      /sciclone/pscr/zxiao01/lib/lapack-3.8.0/liblapack.a \
      /sciclone/pscr/zxiao01/lib/lapack-3.8.0/librefblas.a\
      /sciclone/pscr/zxiao01/lib/gmp-6.1.2/.libs/libgmp.a\
      /sciclone/pscr/zxiao01/lib/lapack-3.8.0/libtmglib.a

FFLAGS= -DPOINTER_SIZE=8 -D MPI

HEADS= -I/sciclone/pscr/zxiao01/lib/sprng2.0/include

MODSF90= sprng_rndm.F90

MODSO=   sprng_rndm.o

# without MODSO
OBJ =   timing_module.o io_module.o matrixcal.o caldet.o matrix_cpmc.o\
        params.o initial_end.o lattice_inform.o initial_kv.o initial_phiT.o matrix_dc.o\
        initial_phi.o modified_GS.o population.o adjustET.o update.o backup_phi.o\
        meas_sisj.o meas_ninj.o meas_pairing.o meas_cicj.o meas_htwo.o measure.o  data_manip.o file_name.o back_propagation.o\
        cpmc.o rcpmc.o step_qmc.o main.o measure_dynamics.o dynamics.o print_observables.o measure_energy.o\
        response_functions.o green_particles.o green_holes.o bcs_update.o bcs_green_pure.o HF.o

#out MODSF90
FUNS=   timing_module.F90 sprng_rndm.F90 io_module.F90 matrixcal.F90 caldet.F90 matrix_cpmc.F90\
        params.F90 initial_end.F90 lattice_inform.F90 initial_kv.F90 initial_phiT.F90 matrix_dc.F90\
        initial_phi.F90 modified_GS.F90 population.F90 adjustET.F90 update.F90 backup_phi.F90\
        meas_sisj.F90 meas_ninj.F90 meas_pairing.F90 meas_cicj.F90 meas_htwo.F90 measure.F90  data_manip.F90 file_name.F90 back_propagation.F90\
        cpmc.F90 rcpmc.F90 step_qmc.F90 main.F90 measure_dynamics.F90 dynamics.F90 print_observables.F90 measure_energy.F90\
        response_functions.F90 green_particles.F90 green_holes.F90 bcs_update.F90 bcs_green_pure.F90 HF.F90


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
	rm -f $(PROG) *.mod *.dat *.out *.o *.info







