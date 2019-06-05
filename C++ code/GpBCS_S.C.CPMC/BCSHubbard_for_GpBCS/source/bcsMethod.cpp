//
// Created by Hao Shi on 12/18/18.
//

#include "../include/bcsMethod.h"

using namespace std;

BcsMethod::BcsMethod()
{

}

BcsMethod::~BcsMethod()
{

}

void BcsMethod::read(const std::string &filename)
{
    ifstream file;
    file.open(filename, ios::in);
    if ( ! file.is_open() ) {cout << "Error opening file in File!!! "<<filename<<endl; exit(1);}

    file>>initialType;
    file>>convergeTolerance;
    file>>maxIterateStep;
    file>>annealMagnitude;
    file>>annealStep;
    file>>relaxMagnitude;
    file>>initMu;
    file>>deltaMu;
    file>>Neta;
    file>>seed;

    file.close();

    analysis();
}

#ifdef MPI_HAO
void MPIBcast(BcsMethod &buffer, int root, MPI_Comm const &comm)
{
    MPIBcast(buffer.initialType, root, comm);
    MPIBcast(buffer.convergeTolerance, root, comm);
    MPIBcast(buffer.maxIterateStep, root, comm);
    MPIBcast(buffer.annealMagnitude, root, comm);
    MPIBcast(buffer.annealStep, root, comm);
    MPIBcast(buffer.relaxMagnitude, root, comm);
    MPIBcast(buffer.initMu, root, comm);
    MPIBcast(buffer.deltaMu, root, comm);
    MPIBcast(buffer.Neta, root, comm);
    MPIBcast(buffer.seed, root, comm);
}
#endif

void BcsMethod::analysis()
{
    if( convergeTolerance < 0.0 )
    {
        cout<<"Errorr!!! convergeTolerance can not be negative!"<<endl;
        exit(1);
    }

    if( relaxMagnitude <= 0.0 )
    {
        cout<<"Errorr!!! relaxMagnitude must be positive!"<<endl;
        exit(1);
    }
}