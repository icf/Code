//
// Created by boruoshihao on 12/25/16.
// Modefied by icf on 04/29/2019
//

#ifndef AFQMCLAB_BCS_H
#define AFQMCLAB_BCS_H

#include "afqmclab.h"

//general BCS

#ifdef MPI_HAO
class BCS;
void MPIBcast(BCS &buffer, int root=0,  const MPI_Comm& comm=MPI_COMM_WORLD);
#endif

class BCS
{
 private:
    std::complex<double> logw;
    tensor_hao::TensorHao<std::complex<double>,2> wf;

 public:
    tensor_hao::TensorHao<std::complex<double>,2> orbital;
    tensor_hao::TensorHao<std::complex<double>,2> occupancy;
//
    BCS();
    BCS(size_t L);
    BCS(const BCS& x);
    BCS(BCS&& x);
    ~BCS();

    BCS & operator  = (const BCS& x);
    BCS & operator  = (BCS&& x);

    const std::complex<double> &getLogw() const;
    const tensor_hao::TensorHao<std::complex<double>, 2> &getWf() const;
    std::complex<double> &logwRef();
    tensor_hao::TensorHao<std::complex<double>, 2> &wfRef();
    size_t getL() const;
//    size_t getN() const;             //icf: without getN

    void resize(size_t L);             //icf: only L=2*Nsite
    void stabilize();
    std::complex<double> normalize();
//    std::complex<double> normalize(double &ratio);
    void addLogw(std::complex<double> logw_add);
    void randomFill();

    void read(const std::string& filename);
    void readOrbital(const std::string& filename);
    void readOccupancy(const std::string& filename);
    void write(const std::string& filename) const;
    int returnNbuf() const;
    double getMemory() const;

#ifdef MPI_HAO
    friend void MPIBcast(BCS &buffer, int root,  const MPI_Comm& comm);
    void pack( std::vector<char> &buf,  int &posit ) const;
    void unpack( const std::vector<char> &buf, int &posit );
#endif

 private:
    void copy_deep(const BCS &x);
    void move_deep(BCS &x);
};

#endif //AFQMCLAB_BCS_H
