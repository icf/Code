//
// Created by boruoshihao on 12/18/8.
//

#include "../include/bcs.h"

using namespace std;
using namespace tensor_hao;

int main(int argc, char** argv)
{
    MPIInit(argc,argv);

    TimerHao timer;
    timer.start(); if(MPIRank()==0) {cout<<"\nStart time:\n";timer.printCurrentTime();}

#ifdef USE_MAGMA
    magma_init();
#endif

    Bcs bcs;
    bcs.run();

#ifdef USE_MAGMA
    magma_finalize();
#endif

    timer.end(); if(MPIRank()==0) {cout<<"\n\nEnd time:\n";timer.printCurrentTime();}
    if(MPIRank()==0) {cout<<"The program's running walltime is:\n"; timer.printFormat();}

    MPIFinalize();
    return 0;
}