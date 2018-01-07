/***************************************************************************/
/*            ____Demonstrates the use of sprng in C++____                 */
/* A random number stream is initialized and a few random double precision */
/* numbers are printed.                                                    */
/***************************************************************************/
#include <stdio.h>
#include <string.h>
#include <iostream.h>

#include "sprng.h"  /* SPRNG header file                                   */

#define SEED 985456376



main()
{
  int streamnum, nstreams, *stream;
  double rn;
  int irn;
  int i;
  int gtype; /*---*/
/*--- reading in a generator type */
#include "gen_types_menu.h"
  printf("Type in a generator type (integers: 0,1,2,3,4,5):  ");
  scanf("%d", &gtype);       
/*---*/  
  /****************** Initialization values *******************************/
            
  streamnum = 0;
  nstreams = 1;

  stream = init_sprng(gtype,streamnum,nstreams,SEED,SPRNG_DEFAULT); /* initialize stream */
  cout << " Print information about new stream:\n";
  print_sprng(stream);	

  /*********************** print random numbers ***************************/

  cout << " Printing 3 random numbers in [0,1):\n";
  for (i=0;i<3;i++)
  {
    rn = sprng(stream);		/* generate a double precision random number */
    cout << rn << "\n";
  }

  free_sprng(stream);  /* free memory used to store stream state */

}
