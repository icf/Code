//
// Created by boruoshihao on 1/13/17.
// Modefied by icf on 04/29/2019
//

#include "../include/HubbardSOCMeasureObserveBCSSD.h"
#include "afqmclab.h"

using namespace std;
using namespace tensor_hao;

HubbardSOCMeasureObserveBCSSD::HubbardSOCMeasureObserveBCSSD()
{
    initModelNullptr();
    reSet();
}

HubbardSOCMeasureObserveBCSSD::HubbardSOCMeasureObserveBCSSD(const HubbardSOC &hubbardSOC_)
{
    setModel( hubbardSOC_ );
    reSet();
}

HubbardSOCMeasureObserveBCSSD::~HubbardSOCMeasureObserveBCSSD()
{

}

void HubbardSOCMeasureObserveBCSSD::reSet()
{
    HubbardSOCMeasureCommuteBCSSD::reSet();

    complex<double> zero(0,0);
    greenMatrixNum = zero;
    densityDensityMatrixNum = zero;
    splusSminusNum = zero;
    sminusSplusNum = zero;
    spairSpairNum = zero;
}

void HubbardSOCMeasureObserveBCSSD::addMeasurement(BCSSDOperation &bcssdOperation, complex<double> denIncrement)
{
    HubbardSOCMeasureCommuteBCSSD::addMeasurement(bcssdOperation, denIncrement);
    addGreenMatrix(bcssdOperation, denIncrement);
    addDensityDensityMatrix(bcssdOperation, denIncrement);
    //addSplusSminus(bcssdOperation, denIncrement);
    //addSminusSplus(bcssdOperation, denIncrement);
    //addSpairSpair(bcssdOperation, denIncrement);
}

void HubbardSOCMeasureObserveBCSSD::addMeasurementBCSBP(BCSSDOperation &bcssdOperation, complex<double> denIncrement)
{
    HubbardSOCMeasureCommuteBCSSD::addMeasurementBCSBP(bcssdOperation, denIncrement);

    addGreenMatrixBCSBP(bcssdOperation, denIncrement);
    addDensityDensityMatrixBCSBP(bcssdOperation, denIncrement);
    //addSplusSminus(bcssdOperation, denIncrement);
    //addSminusSplus(bcssdOperation, denIncrement);
    //addSpairSpair(bcssdOperation, denIncrement);

}

void HubbardSOCMeasureObserveBCSSD::write() const
{
    HubbardSOCMeasureCommuteBCSSD::write();
    writeKNumVumRum();
    writeThreadSum(greenMatrixNum.size(), greenMatrixNum.data(), "greenMatrixNum.dat", ios::app);
    writeThreadSum(densityDensityMatrixNum.size(), densityDensityMatrixNum.data(), "densityDensityMatrixNum.dat", ios::app);
    //writeThreadSum(splusSminusNum.size(), splusSminusNum.data(), "splusSminusNum.dat", ios::app);
    //writeThreadSum(sminusSplusNum.size(), sminusSplusNum.data(), "sminusSplusNum.dat", ios::app);
    //writeThreadSum(spairSpairNum.size(), spairSpairNum.data(), "spairSpairNum.dat", ios::app);
}

void HubbardSOCMeasureObserveBCSSD::mixedWrite() const
{
    HubbardSOCMeasureCommuteBCSSD::mixedWrite();
}

double HubbardSOCMeasureObserveBCSSD::getMemory() const
{
    double mem(0.0);
    mem += HubbardSOCMeasureCommuteBCSSD::getMemory();
    mem += greenMatrixNum.getMemory();
    mem += densityDensityMatrixNum.getMemory();
    mem += splusSminusNum.getMemory();
    mem += sminusSplusNum.getMemory();
    mem += spairSpairNum.getMemory();

    return mem;
}

void HubbardSOCMeasureObserveBCSSD::addGreenMatrix(BCSSDOperation &bcssdOperation, complex<double> denIncrement)
{
    const TensorHao< complex<double>, 2 > &greenMatrix = bcssdOperation.returnGreenMatrix();

    size_t L2 = getHubbardSOC()->getL() * 2;

    if( greenMatrixNum.rank(0) != L2 ) { greenMatrixNum.resize(L2, L2); greenMatrixNum = complex<double>(0,0); }

    greenMatrixNum += ( greenMatrix * denIncrement );
}

void HubbardSOCMeasureObserveBCSSD::addGreenMatrixBCSBP(BCSSDOperation &bcssdOperation, complex<double> denIncrement)
{
    const TensorHao< complex<double>, 2 > &greenMatrix = bcssdOperation.returnGreenMatrixBCSBP();

    size_t L2 = getHubbardSOC()->getL() * 2;

    if( greenMatrixNum.rank(0) != L2 ) { greenMatrixNum.resize(L2, L2); greenMatrixNum = complex<double>(0,0); }

    greenMatrixNum += ( greenMatrix * denIncrement );
}

void HubbardSOCMeasureObserveBCSSD::addDensityDensityMatrix(BCSSDOperation &bcssdOperation, complex<double> denIncrement)
{
    const TensorHao< complex<double>, 2 > &densityDensityMatrix = bcssdOperation.returnDensityDensityMatrix();

    size_t L2 = getHubbardSOC()->getL() * 2;

    if( densityDensityMatrixNum.rank(0) != L2 ) { densityDensityMatrixNum.resize(L2, L2); densityDensityMatrixNum = complex<double>(0,0); }

    densityDensityMatrixNum += ( densityDensityMatrix * denIncrement );
}

void HubbardSOCMeasureObserveBCSSD::addDensityDensityMatrixBCSBP(BCSSDOperation &bcssdOperation, complex<double> denIncrement)
{
    const TensorHao< complex<double>, 2 > &densityDensityMatrix = bcssdOperation.returnDensityDensityMatrixBCSBP();

    size_t L2 = getHubbardSOC()->getL() * 2;

    if( densityDensityMatrixNum.rank(0) != L2 ) { densityDensityMatrixNum.resize(L2, L2); densityDensityMatrixNum = complex<double>(0,0); }

    densityDensityMatrixNum += ( densityDensityMatrix * denIncrement );
}

const TensorHao<complex<double>, 2> HubbardSOCMeasureObserveBCSSD::returnGreenMatrix()
{
    size_t L2 = getHubbardSOC()->getL() * 2;

    TensorHao<complex<double>, 2> greenMatrix;
    TensorHao<complex<double>, 2> greenMatrixNumTot(L2,L2);
    MPISum(L2*L2,greenMatrixNum.data(),greenMatrixNumTot.data());

    complex<double> denTot = MPISum(den);

    greenMatrix = greenMatrixNumTot/denTot;

    //MPIBcast(greenMatrix);  since we don't do Bcast here, only greenMatrix at Rank == 0 is valid

    return greenMatrix;
}

const TensorHao<complex<double>, 2> HubbardSOCMeasureObserveBCSSD::returnDensityDensityMatrix()
{
    size_t L2 = getHubbardSOC()->getL() * 2;

    TensorHao<complex<double>, 2> densityDensityMatrix;
    TensorHao<complex<double>, 2> densityDensityMatrixNumTot(L2,L2);
    MPISum(L2*L2,densityDensityMatrixNum.data(),densityDensityMatrixNumTot.data());

    complex<double> denTot = MPISum(den);

    densityDensityMatrix = densityDensityMatrixNumTot/denTot;
    //MPIBcast(greenMatrix);  since we don't do Bcast here, only greenMatrix at Rank == 0 is valid

    return densityDensityMatrix;
}

HubbardSOCMeasureObserveBCSSD::HubbardSOCMeasureObserveBCSSD(const HubbardSOCMeasureObserveBCSSD &x)
{

}

HubbardSOCMeasureObserveBCSSD &HubbardSOCMeasureObserveBCSSD::operator=(const HubbardSOCMeasureObserveBCSSD &x)
{
    return *this;
}
