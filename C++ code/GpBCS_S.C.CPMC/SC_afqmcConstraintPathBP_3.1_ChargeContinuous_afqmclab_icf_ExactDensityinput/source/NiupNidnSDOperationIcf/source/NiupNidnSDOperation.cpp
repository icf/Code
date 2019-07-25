//
// Created by boruoshihao on 1/10/17.
//

#include "../include/NiupNidnSDOperationIcf.h"
#include "afqmclab.h"

using namespace std;
using namespace tensor_hao;

NiupNidnSampleSDOperationIcf::NiupNidnSampleSDOperationIcf()  { }

NiupNidnSampleSDOperationIcf::~NiupNidnSampleSDOperationIcf() { }

void NiupNidnSampleSDOperationIcf::applyToRight(const NiupNidnSampleIcf &oneBody, const SD &walker, SD &walkerNew) const
{
    size_t L = walker.getL(); size_t N = walker.getN(); size_t halfL = oneBody.getL();

    if( L != halfL*2 ) { cout<<"Error!!! NiupNidnSampleIcf size is not consistent with walker!"<<endl; exit(1); }
    if( walkerNew.getL() != L  ||  walkerNew.getN() != N ) walkerNew.wfRef().resize( L, N );

    const TensorHao<complex<double>,2> &wf = walker.getWf();
    TensorHao<complex<double>,2> &wfNew = walkerNew.wfRef();

    const TensorHao<complex<double>,1> &diag00 = oneBody.diag00;
    const TensorHao<complex<double>,1> &diag10 = oneBody.diag10;
    const TensorHao<complex<double>,1> &diag01 = oneBody.diag01;
    const TensorHao<complex<double>,1> &diag11 = oneBody.diag11;

    for(size_t j = 0; j < N; ++j)
    {
        for(size_t i = 0; i < halfL; ++i)
        {
            wfNew(i,j)        = diag00(i) * wf(i,j) + diag01(i)*wf(i+halfL, j);
            wfNew(i+halfL, j) = diag10(i) * wf(i,j) + diag11(i)*wf(i+halfL, j);
        }
    }

    walkerNew.logwRef() = oneBody.logw + walker.getLogw();

}

void NiupNidnSampleSDOperationIcf::applyToLeft(const NiupNidnSampleIcf &oneBody, const SD &walker, SD &walkerNew) const
{
    size_t L = walker.getL(); size_t N = walker.getN(); size_t halfL = oneBody.getL();

    if( L != halfL*2 ) { cout<<"Error!!! NiupNidnSampleIcf size is not consistent with walker!"<<endl; exit(1); }
    if( walkerNew.getL() != L  ||  walkerNew.getN() != N ) walkerNew.wfRef().resize( L, N );

    const TensorHao<complex<double>,2> &wf = walker.getWf();
    TensorHao<complex<double>,2> &wfNew = walkerNew.wfRef();

    const TensorHao<complex<double>,1> diag00 = conj( oneBody.diag00 );
    const TensorHao<complex<double>,1> diag10 = conj( oneBody.diag01 );
    const TensorHao<complex<double>,1> diag01 = conj( oneBody.diag10 );
    const TensorHao<complex<double>,1> diag11 = conj( oneBody.diag11 );

    for(size_t j = 0; j < N; ++j)
    {
        for(size_t i = 0; i < halfL; ++i)
        {
            wfNew(i,j)        = diag00(i) * wf(i,j) + diag01(i)*wf(i+halfL, j);
            wfNew(i+halfL, j) = diag10(i) * wf(i,j) + diag11(i)*wf(i+halfL, j);
        }
    }

    walkerNew.logwRef() = conj( oneBody.logw ) + walker.getLogw();
}

