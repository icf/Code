

       program simplef_simple
       implicit none

#define SIMPLE_SPRNG	
#include "sprng_f.h"

       integer i
       real*8 rn
       
       print *, 'Printing 3 random numbers in [0,1):'
       do 100 i = 1, 3
          rn = sprng()
          write(*,"(f8.6)")  rn
 100   continue


       end
