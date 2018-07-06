FC = ifort
#FC=pgf90
#LIBS = 
PROG = tcpqmc
#--------- SUN 
LIBS = -L/pub/opt/lib -llapack -lblas 
#LIBS     = -L ../usr/local/lib -xlic_lib=imsl
#FFLAGS = -dalign -fast
#FFLAGS =  -Mmpi
#FFLAGS = -openmp

#FFLAGS =  -g -traceback -check all # -warn all  #-O3 -ffree-line-length-none #Set the optimization flags

#FFLAGS =  -O3 #-check all #-warn all  #-O3 -ffree-line-length-none #Set the optimization flags
#--------- DEC
#FFLAGS = -g
#FFLAGS = -migrate -tune host -O5 -fast

OBJ = params.o io_module.o matrixcal.o caldet.o zqr.o \
      params1.o Hopt_module.o calphy_module.o useva05_2.o va05.o diaMat.o\
      gasdev.o rannyu.o  \
      initialize.o ini_EDdynamics.o\
      cpqmc_subs.o\
      reconfiguration.o \
      forward_propagation.o \
      lanczos_onestep.o sub_lanc.o hubbard_chain_subs.o buildbasis.o\
      slatertofock.o perm_cy_unordered.o \
      main_forward.o \
      ED_onestep.o assert.f90

     


all: $(PROG)

clobber:
	rm -f *.o *.M *.mod

clean:clobber
	rm -f mf1

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

