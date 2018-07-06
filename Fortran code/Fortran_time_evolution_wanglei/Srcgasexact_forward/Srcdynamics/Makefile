FC = ifort
#FC=pgf90
#LIBS = 
PROG = hubbardchain
#--------- SUN 
LIBS = -L/pub/opt/lib -llapack -lblas 

FFLAGS =  #-check all #-warn all  #-O3 -ffree-line-length-none #Set the optimization flags
#LIBS     = -L ../usr/local/lib -xlic_lib=imsl
#FFLAGS = -dalign -fast
#FFLAGS =  -Mmpi
#FFLAGS = -openmp
#--------- DEC
#FFLAGS = -g
#FFLAGS = -migrate -tune host -O5 -fast

OBJ = matrixcal.o params.o calphy_module.o Hopt_module.o io_module.o\
      buildbasis.o  hubbard_chain_subs.o \
      diaMat.o sub_lanc.o \
      slatertofock.o perm_cy_unordered.o caldet.o \
      lanczos_onestep.f90\
      readindata.o \


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

