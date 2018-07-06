FC = ifort
#FC=pgf90
#LIBS = 
PROG = slaterfit
#--------- SUN 
LIBS = -L/pub/opt/lib -llapack -lblas 

FFLAGS =  -g -traceback -check all #-warn all  #-O3 -ffree-line-length-none #Set the optimization flags
#LIBS     = -L ../usr/local/lib -xlic_lib=imsl
#FFLAGS = -dalign -fast
#FFLAGS =  -Mmpi
#FFLAGS = -openmp
#--------- DEC
#FFLAGS = -g
#FFLAGS = -migrate -tune host -O5 -fast

OBJ = useva05.o va05.o rannyu.o caldet.o matrixcal.o main.o 


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

