FC = g++

PROG = VMpBCS

LIBS= 

FFLAGS= 

HEADS= -I/home/xiaozhiyu/Documents/LBFGSpp-master/include \
       -I/home/xiaozhiyu/Documents/eigen

OBJ =   VMpBCS.o 

FUNS=   VMpBCS.cpp 


$(PROG): $(OBJ) 
	$(FC) $(FFLAGS) $(HEADS) -o $(PROG) $(OBJ) $(LIBS)

$(OBJ): 
	$(FC) -c $(FFLAGS) $(HEADS) $(FUNS) $(LIBS)

	

.PHONY : dataclean refrash realclean

dateclean:
	rm -f *.dat

refrash:
	rm -f *.mod
	$(PROG)

realclean:
	rm -f $(PROG) *.mod *.dat *.out *.o *.info







