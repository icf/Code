//
// Created by boruoshihao on 1/10/17.
// Modefied by icf on 05/06/2019
//

#include "../include/NiupNidnBCSOperation.h"
#include "afqmclab.h"

using namespace std;
using namespace tensor_hao;

NiupNidnSampleBCSOperation::NiupNidnSampleBCSOperation()  { }

NiupNidnSampleBCSOperation::~NiupNidnSampleBCSOperation() { }

void NiupNidnSampleBCSOperation::applyToRight(const NiupNidnSample &oneBody, const BCS &walker, BCS &walkerNew) const
{
    size_t L = walker.getL(); size_t halfL = oneBody.getL();

    if( L != halfL*2 ) { cout<<"Error!!! NiupNidnSample size is not consistent with walker!"<<endl; exit(1); }
    if( walkerNew.getL() != L  ) walkerNew.wfRef().resize( L,L );

    const TensorHao<complex<double>,2> &wf = walker.getWf();
    TensorHao<complex<double>,2> &wfNew = walkerNew.wfRef();\
    TensorHao< complex<double>, 2 > tempOne(L,L);tempOne=complex <double> (0,0);

    const TensorHao<complex<double>,1> &diag00 = oneBody.diag00;
    const TensorHao<complex<double>,1> &diag10 = oneBody.diag10;
    const TensorHao<complex<double>,1> &diag01 = oneBody.diag01;
    const TensorHao<complex<double>,1> &diag11 = oneBody.diag11;

    for(size_t j = 0; j < L; ++j)
    {
        for(size_t i = 0; i < halfL; ++i)
        {
            tempOne(i,j)        = diag00(i) * wf(i,j) + diag01(i)*wf(i+halfL, j);
            tempOne(i+halfL, j) = diag10(i) * wf(i,j) + diag11(i)*wf(i+halfL, j);
        }
    }
    for(size_t j = 0; j < halfL; ++j)
    {
        for(size_t i = 0; i < L; ++i)
        {
            wfNew(i,j)        = diag00(j) * tempOne(i,j) + diag10(j)*tempOne(i, j+halfL);
            wfNew(i,j+halfL)  = diag01(j) * tempOne(i,j) + diag11(j)*tempOne(i, j+halfL);
        }
    }

    for(size_t j = 0; j < L; ++j)
    {
        for(size_t i = 0; i < halfL; ++i)
        {
            walkerNew.orbital(i,j)        = diag00(i) * walker.orbital(i,j) + diag01(i)*walker.orbital(i+halfL, j);
            walkerNew.orbital(i+halfL, j) = diag10(i) * walker.orbital(i,j) + diag11(i)*walker.orbital(i+halfL, j);
        }
    }
    walkerNew.occupancy=walker.occupancy;
    walkerNew.logwRef() = oneBody.logw + walker.getLogw();

}

void NiupNidnSampleBCSOperation::applyToLeft(const NiupNidnSample &oneBody, const BCS &walker, BCS &walkerNew) const
{
    size_t L = walker.getL(); size_t halfL = oneBody.getL();

    if( L != halfL*2 ) { cout<<"Error!!! NiupNidnSample size is not consistent with walker!"<<endl; exit(1); }

    if( walkerNew.getL() != L ) walkerNew.resize( L );

    const TensorHao<complex<double>,2> &wf = walker.getWf();
    TensorHao<complex<double>,2> &wfNew = walkerNew.wfRef();
    TensorHao< complex<double>, 2 > tempOne(L,L);tempOne=complex <double> (0,0);

    const TensorHao<complex<double>,1> diag00 = conj( oneBody.diag00 );
    const TensorHao<complex<double>,1> diag10 = conj( oneBody.diag01 );
    const TensorHao<complex<double>,1> diag01 = conj( oneBody.diag10 );
    const TensorHao<complex<double>,1> diag11 = conj( oneBody.diag11 );

    for(size_t j = 0; j < L; ++j)
    {
        for(size_t i = 0; i < halfL; ++i)
        {
            tempOne(i,j)        = diag00(i) * wf(i,j) + diag01(i)*wf(i+halfL, j);
            tempOne(i+halfL, j) = diag10(i) * wf(i,j) + diag11(i)*wf(i+halfL, j);
        }
    }
    for(size_t j = 0; j < halfL; ++j)
    {
        for(size_t i = 0; i < L; ++i)
        {
            wfNew(i,j)        = diag00(j) * tempOne(i,j) + diag10(j)*tempOne(i, j+halfL);
            wfNew(i,j+halfL)  = diag01(j) * tempOne(i,j) + diag11(j)*tempOne(i, j+halfL);
        }
    }

    for(size_t j = 0; j < L; ++j)
    {
        for(size_t i = 0; i < halfL; ++i)
        {
            walkerNew.orbital(i,j)        = diag00(i) * walker.orbital(i,j) + diag01(i)*walker.orbital(i+halfL, j);
            walkerNew.orbital(i+halfL, j) = diag10(i) * walker.orbital(i,j) + diag11(i)*walker.orbital(i+halfL, j);
        }
    }  
    walkerNew.occupancy=walker.occupancy;
    walkerNew.logwRef() = conj( oneBody.logw ) + walker.getLogw();
}

