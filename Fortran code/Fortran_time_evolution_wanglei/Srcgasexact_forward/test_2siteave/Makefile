FC = ifort
#FC=pgf90
#LIBS = 
PROG = 2site
#--------- SUN 
LIBS = -L/pub/opt/lib -llapack -lblas 
#LIBS     = -L ../usr/local/lib -xlic_lib=imsl
#FFLAGS = -dalign -fast
#FFLAGS =  -Mmpi
#FFLAGS = -openmp

#FFLAGS =  -g -traceback #-check all# -warn all  #-O3 -ffree-line-length-none #Set the optimization flags

#FFLAGS =  -O3 #-check all #-warn all  #-O3 -ffree-line-length-none #Set the optimization flags
#--------- DEC
#FFLAGS = -g
#FFLAGS = -migrate -tune host -O5 -fast

OBJ = matrixcal.o 2siteed_ave.o     


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

