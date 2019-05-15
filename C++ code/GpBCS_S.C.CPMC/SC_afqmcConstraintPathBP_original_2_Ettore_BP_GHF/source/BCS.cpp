//
// Created by boruoshihao on 12/25/16.
// Modefied by icf on 04/29/2019
//

#include "../include/BCS.h"
#include "afqmclab.h"

using namespace std;
using namespace tensor_hao;

BCS::BCS():logw(0.0) {}

BCS::BCS(size_t L) :logw(0.0) { wf.resize(L, L);orbital.resize(L, L);occupancy.resize(L, L); }    //icf: general BCS is 2*Nsite by 2*Nsite matrix

BCS::BCS(const BCS &x) { copy_deep(x); }

BCS::BCS(BCS &&x) { move_deep(x); }

BCS::~BCS() { }

BCS &BCS::operator=(const BCS &x) { copy_deep(x); return *this; }

BCS &BCS::operator=(BCS &&x) { move_deep(x); return *this; }

const complex<double> &BCS::getLogw() const { return logw; }

const TensorHao<complex<double>, 2> &BCS::getWf() const { return wf; }

complex<double> &BCS::logwRef() { return logw; }

TensorHao<complex<double>, 2> &BCS::wfRef() { return wf; }  //icf: don't understand the logic of this.????????????? why not just give public wf?

size_t BCS::getL() const { return wf.rank(0); }

void BCS::resize(size_t L)
{
    wf.resize(L, L);
    orbital.resize(L, L);
    occupancy.resize(L, L);
}

void BCS::stabilize()
{
    size_t L=getL();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    TensorHao<double, 1> orbitalR(L);
    TensorHao<complex<double>, 2> orbitalRMatrix(L,L);
    BL_NAME(QRMatrix)(orbital,orbitalR);
    for(size_t i=1-1;i<=L-1;i++){
       orbitalRMatrix(i,i)=complex <double>(orbitalR(i));
    }

    TensorHao<complex<double>, 2> tempMatrixOne(L,L);
    TensorHao<complex<double>, 2> tempMatrixOneOne(L,L);
    BL_NAME(gmm)(orbitalRMatrix, occupancy, tempMatrixOne ); 
    BL_NAME(gmm)(tempMatrixOne, orbitalRMatrix, tempMatrixOneOne ,'N','T');
    occupancy=tempMatrixOneOne;

    complex <double> sumTemp;
    for(size_t i=1-1;i<=L-1;i++){
    for(size_t j=1-1;j<=L-1;j++){
       sumTemp += abs(occupancy(i,j));
    }
    }
    occupancy = occupancy/sumTemp;
    logw += log(sumTemp);
    TensorHao<complex<double>, 2> tempMatrixWf(L,L);
    TensorHao<complex<double>, 2> tempMatrixTwo(L,L);
    BL_NAME(gmm)(orbital, occupancy, tempMatrixTwo ); 
    BL_NAME(gmm)(tempMatrixTwo, orbital, tempMatrixWf,'N','T'); 
    wf=tempMatrixWf;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
}

complex<double> BCS::normalize()
{
    stabilize();
    complex<double> logwTemp(logw);
    logw=0.0;
    return logwTemp;
}

//std::complex<double> BCS::normalize(double &ratio)
//{
//    complex<double> logwTemp(logw);
//    logw=0.0;
//    return logwTemp;
//}

void BCS::addLogw(std::complex<double> logw_add)
{
    logw += logw_add;
}

void BCS::randomFill()
{
    //tensor_hao::randomFill(wf);
    tensor_hao::randomFill(orbital);
    tensor_hao::randomFill(occupancy);
    normalize();
}

void BCS::read(const string &filename)    //icf: Be careful about the read and write form
{
    ifstream file;
    file.open(filename, ios::in);
    if ( ! file.is_open() ) { cout << "Error opening file in File!!! "<<filename<<endl; exit(1); }
    readFile(logw, file);
    wf.read(file);
    file.close();
}
void BCS::readOrbital(const string &filename)    //icf: Be careful about the read and write form
{
    ifstream file;
    file.open(filename, ios::in);
    if ( ! file.is_open() ) { cout << "Error opening file in File!!! "<<filename<<endl; exit(1); }
    readFile(logw, file);
    orbital.read(file);
    file.close();
}
void BCS::readOccupancy(const string &filename)    //icf: Be careful about the read and write form
{
    ifstream file;
    file.open(filename, ios::in);
    if ( ! file.is_open() ) { cout << "Error opening file in File!!! "<<filename<<endl; exit(1); }
    readFile(logw, file);
    occupancy.read(file);
    file.close();
}

void BCS::write(const string &filename) const
{
    ofstream file;
    file.open(filename, ios::out|ios::trunc);
    if ( ! file.is_open() ) { cout << "Error opening file in File!!! "<<filename<<endl; exit(1); }
    writeFile(logw, file);
    wf.write(file);
    file.close();
}

int BCS::returnNbuf() const
{
    return 16+16*wf.size()+16*orbital.size()+16*occupancy.size();
}

double BCS::getMemory() const
{
    return 16.0+wf.getMemory()+16*orbital.getMemory()+16*occupancy.getMemory();
}

#ifdef MPI_HAO
void MPIBcast(BCS &buffer, int root, MPI_Comm const &comm)
{
    MPIBcast( buffer.logw, root, comm );
    MPIBcast( buffer.wf, root, comm );
    MPIBcast( buffer.orbital, root, comm );
    MPIBcast( buffer.occupancy, root, comm );
}

void BCS::pack(vector<char> &buf, int &posit) const
{
    MPI_Pack(&logw, 1, MPI_DOUBLE_COMPLEX, buf.data(), buf.size(), &posit, MPI_COMM_WORLD);
    MPI_Pack(wf.data(), wf.size(), MPI_DOUBLE_COMPLEX, buf.data(), buf.size(), &posit, MPI_COMM_WORLD);
    MPI_Pack(orbital.data(), orbital.size(), MPI_DOUBLE_COMPLEX, buf.data(), buf.size(), &posit, MPI_COMM_WORLD);
    MPI_Pack(occupancy.data(), occupancy.size(), MPI_DOUBLE_COMPLEX, buf.data(), buf.size(), &posit, MPI_COMM_WORLD);
}

void BCS::unpack(const vector<char> &buf, int &posit)
{
    MPI_Unpack(buf.data(), buf.size(), &posit, &logw, 1, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
    MPI_Unpack(buf.data(), buf.size(), &posit, wf.data(), wf.size(), MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
    MPI_Unpack(buf.data(), buf.size(), &posit, orbital.data(), orbital.size(), MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
    MPI_Unpack(buf.data(), buf.size(), &posit, occupancy.data(), occupancy.size(), MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
}
#endif

void BCS::copy_deep(const BCS &x)
{
    logw = x.logw;
    wf = x.wf;
    orbital = x.orbital;
    occupancy = x.occupancy;
}

void BCS::move_deep(BCS &x)
{
    logw = x.logw;
    wf = move( x.wf );
    orbital = move( x.orbital );
    occupancy = move( x.occupancy );
}


