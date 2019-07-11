//
// Created by boruoshihao on 12/28/16.
// Modefied by icf on 05/06/2019
//

#include "../include/hopBCSOperation.h"

using namespace std;
using namespace tensor_hao;

HopBCSOperation::HopBCSOperation() { }

HopBCSOperation::~HopBCSOperation() { }

void HopBCSOperation::applyToRight(const Hop &oneBody, const BCS &walker, BCS &walkerNew) const
{
    checkAndResize(oneBody, walker, walkerNew);
    size_t L = walker.getL();
    TensorHao< complex<double>, 2 > tempOne(L,L);tempOne=complex <double> (0,0);
    BL_NAME(gmm)( oneBody.matrix, walker.getWf(), tempOne );
    BL_NAME(gmm)( tempOne, oneBody.matrix, walkerNew.wfRef(),'N','T' );

    BL_NAME(gmm)( oneBody.matrix, walker.orbital, walkerNew.orbital );

    walkerNew.occupancy=walker.occupancy;
    walkerNew.logwRef() = oneBody.logw + walker.getLogw();
}

void HopBCSOperation::applyToLeft(const Hop &oneBody, const BCS &walker, BCS &walkerNew) const
{
    checkAndResize(oneBody, walker, walkerNew);
    size_t L = walker.getL();
    TensorHao< complex<double>, 2 > tempOne(L,L);tempOne=complex <double> (0,0);
    BL_NAME(gmm)( oneBody.matrix, walker.getWf(), tempOne,'C' );
    BL_NAME(gmm)( tempOne, conj(oneBody.matrix), walkerNew.wfRef() );

    BL_NAME(gmm)( oneBody.matrix, walker.orbital, walkerNew.orbital,'C'  );

    walkerNew.occupancy=walker.occupancy;
    walkerNew.logwRef() = conj( oneBody.logw ) + walker.getLogw();
}

void HopBCSOperation::checkAndResize(const Hop &oneBody, const BCS &walker, BCS &walkerNew) const
{
    size_t L = walker.getL();
    if( oneBody.getL() !=  L ) {cout<<"Error!!! Hop size is not consistent with walker!"<<endl; exit(1); }
    if( walkerNew.getL() != L  ) walkerNew.resize( L );

}
