//
// Created by Hao Shi on 12/18/18.
//

#ifndef BCSHUBBARD_BCSMETHOD_H
#define BCSHUBBARD_BCSMETHOD_H

#include "afqmclab.h"

class BcsMethod
{
 public:
    std::string initialType; //"setFromModel", "readWaveFunction", "readOrderParameter"
    double convergeTolerance;
    size_t maxIterateStep;
    double annealMagnitude;
    size_t annealStep;
    double relaxMagnitude; // 1.0 fully relax to new order paramter, 0.0 not update
    double initMu;
    double deltaMu;
    double Neta;
    int seed;  // -1. read file, 0. random, else is seeds

    BcsMethod();
    ~BcsMethod();

    void read(const std::string& filename);

#ifdef MPI_HAO
    friend void MPIBcast(BcsMethod &buffer, int root=0,  const MPI_Comm& comm=MPI_COMM_WORLD);
#endif

 private:
    void analysis();
};

#endif //BCSHUBBARD_BCSMETHOD_H
