//
// Created by boruoshihao on 1/13/17.
// Modefied by icf on 04/29/2019
// Modefied by icf on 06/13/2019
//

#include "../include/HubbardSOCMeasureCommuteBCSSD.h"
//#include "../../../utilities/manipulateMCData/include/writeThreadSum.h"
#include "afqmclab.h"

using namespace std;
using namespace tensor_hao;

HubbardSOCMeasureCommuteBCSSD::HubbardSOCMeasureCommuteBCSSD()
{
    initModelNullptr();
    reSet();
}

HubbardSOCMeasureCommuteBCSSD::HubbardSOCMeasureCommuteBCSSD(const HubbardSOC &hubbardSOC_)
{
    setModel(hubbardSOC_);
    reSet();
}

HubbardSOCMeasureCommuteBCSSD::~HubbardSOCMeasureCommuteBCSSD()
{

}

const HubbardSOC *HubbardSOCMeasureCommuteBCSSD::getHubbardSOC() const
{
    return hubbardSOC;
}

void HubbardSOCMeasureCommuteBCSSD::initModelNullptr()
{
    hubbardSOC = nullptr;
}

void HubbardSOCMeasureCommuteBCSSD::setModel(const HubbardSOC &hubbardSOC_)
{
    hubbardSOC = &hubbardSOC_;
}

void HubbardSOCMeasureCommuteBCSSD::reSet()
{
    complex<double> zero(0,0);
    den = zero;
    HNum = zero;
    KNum = zero;
    VNum = zero;
    RNum = zero;
}

complex<double> HubbardSOCMeasureCommuteBCSSD::returnEnergy()
{
    complex<double> Vtot   = MPISum(VNum);
    complex<double> Ktot   = MPISum(KNum);
    complex<double> Htot   = MPISum(HNum);
    complex<double> denTot = MPISum(den);
    complex<double> energy;
    if( MPIRank() == 0 ) energy = Htot/denTot;
    if( MPIRank() == 0 ){
       cout<<"Energy: "<<energy<<endl;
       cout<<"K: "<<Ktot/denTot<<endl;
       cout<<"V: "<<Vtot/denTot<<endl;
       cout<<"denTot: "<<denTot<<endl;
    }
    MPIBcast(energy);
    return energy;
}

void HubbardSOCMeasureCommuteBCSSD::addMeasurement(BCSSDOperation &bcssdOperation, complex<double> denIncrement)
{
    den += denIncrement;

    const TensorHao< complex<double>, 2 > &greenMatrix = bcssdOperation.returnGreenMatrix();
    const TensorHao< complex<double>, 2 > &cmcmMatrix = bcssdOperation.returnCmCmMatrix();
    const TensorHao< complex<double>, 2 > &cpcpMatrix = bcssdOperation.returnCpCpMatrix();

    addEnergy(greenMatrix, cmcmMatrix, cpcpMatrix, denIncrement);
}

void HubbardSOCMeasureCommuteBCSSD::addMeasurementBCSBP(BCSSDOperation &bcssdOperation, complex<double> denIncrement)
{
    den += denIncrement;

    const TensorHao< complex<double>, 2 > &greenMatrix = bcssdOperation.returnGreenMatrixBCSBP();
    const TensorHao< complex<double>, 2 > &cmcmMatrix = bcssdOperation.returnCmCmMatrixBCSBP();
    const TensorHao< complex<double>, 2 > &cpcpMatrix = bcssdOperation.returnCpCpMatrixBCSBP();
    
    addEnergy(greenMatrix, cmcmMatrix, cpcpMatrix, denIncrement);

}

NiupNidnForce HubbardSOCMeasureCommuteBCSSD::getForce(const NiupNidn &niupNidn, BCSSDOperation &bcssdOperation, double cap)
{
    size_t halfL = niupNidn.getL(); const string &decompType = niupNidn.getDecompType();

    TensorHao< complex<double>, 1 > backGround(halfL);
    if( decompType == "densityCharge" )
    {
        const TensorHao< complex<double>, 1 > &greenDiagonal = bcssdOperation.returnGreenDiagonal();
        for(size_t i = 0; i < halfL; ++i) backGround(i) = greenDiagonal(i) + greenDiagonal(i+halfL) -1.0;
    }
    else if( decompType == "densitySpin" )
    {
        const TensorHao< complex<double>, 1 > &greenDiagonal = bcssdOperation.returnGreenDiagonal();
        for(size_t i = 0; i < halfL; ++i) backGround(i) = greenDiagonal(i) - greenDiagonal(i+halfL);
    }
    else if( decompType == "hopCharge" )
    {
        const TensorHao< complex<double>, 1 > &greenOffDiagonal = bcssdOperation.returnGreenOffDiagonal();
        for(size_t i = 0; i < halfL; ++i) backGround(i) = greenOffDiagonal(i) + greenOffDiagonal(i+halfL);
    }
    else if( decompType == "hopSpin" )
    {
        const TensorHao< complex<double>, 1 > &greenOffDiagonal = bcssdOperation.returnGreenOffDiagonal();
        for(size_t i = 0; i < halfL; ++i) backGround(i) = greenOffDiagonal(i) - greenOffDiagonal(i+halfL);
    }
    else if( decompType == "chargeDec" )
    {
        complex <double> chargeDecAlpha; chargeDecAlpha=1.2;   //icf: We choose it to be 1 and also set it as 1.0 at NiupNidn.cpp
        const TensorHao< complex<double>, 1 > &greenDiagonal = bcssdOperation.returnGreenDiagonal();
        for(size_t i = 0; i < halfL; ++i) backGround(i) = greenDiagonal(i) + greenDiagonal(i+halfL) - chargeDecAlpha;
    }
    else
    {
        cout<<"Error! Can not find the matched decompType! "<<decompType<<endl;
        exit(1);
    }

    NiupNidnForce force(halfL);

    const TensorHao<complex<double>, 1> &gamma = niupNidn.getGamma();
    if( decompType != "chargeDec" )
    {
       for (size_t i = 0; i < halfL; ++i)
       {
           force(i) = gamma(i) * backGround(i);
           //if( force(i).real() >  cap ) force(i) =  cap;
           //if( force(i).real() < -cap ) force(i) = -cap;
       }
    }else if( decompType == "chargeDec" ){
       for (size_t i = 0; i < halfL; ++i)
       {
           force(i) = gamma(i) * backGround(i);           
           //if( force(i).imag() >  cap ) force(i) = complex <double> (0.0,cap) ;
           //if( force(i).imag() < -cap ) force(i) = complex <double> (0.0,-cap) ;
       }
    }


    return force;
}

void HubbardSOCMeasureCommuteBCSSD::write() const
{
    writeThreadSum(den, "den.dat", ios::app);
    writeThreadSum(HNum, "HNum.dat", ios::app);
}

void HubbardSOCMeasureCommuteBCSSD::mixedWrite() const
{
    writeThreadSum(den, "denMixed.dat", ios::app);
    writeThreadSum(HNum, "HNumMixed.dat", ios::app);
}

void HubbardSOCMeasureCommuteBCSSD::writeKNumVumRum() const
{
    writeThreadSum(KNum, "KNum.dat", ios::app);
    writeThreadSum(VNum, "VNum.dat", ios::app);
    writeThreadSum(RNum, "RNum.dat", ios::app);
}

double HubbardSOCMeasureCommuteBCSSD::getMemory() const
{
    return 8.0+16.0*5;
}

HubbardSOCMeasureCommuteBCSSD::HubbardSOCMeasureCommuteBCSSD(const HubbardSOCMeasureCommuteBCSSD &x)
{

}

HubbardSOCMeasureCommuteBCSSD &HubbardSOCMeasureCommuteBCSSD::operator=(const HubbardSOCMeasureCommuteBCSSD &x)
{
    return *this;
}

void HubbardSOCMeasureCommuteBCSSD::addEnergy(const TensorHao<complex<double>, 2> &greenMatrix, const TensorHao<complex<double>, 2> &cmcmMatrix, const TensorHao<complex<double>, 2> &cpcpMatrix, complex<double> denIncrement)
{
    complex<double> Kenergy(0,0), Venergy(0,0), Renergy(0,0);

    size_t L  = hubbardSOC->getL(); size_t L2 = L*2;

    //Add K
    const TensorHao< complex<double>, 2 > &K = hubbardSOC->getK();
    for(size_t i = 0; i < L2; ++i)
    {
        for(size_t j = 0; j < L2; ++j)
        {
            Kenergy += K(j,i) * greenMatrix(j,i);
        }
    }

    //Add U
    const TensorHao< double, 1> &U = hubbardSOC->getU();
    for(size_t i = 0; i < L; ++i)
    {
        Venergy += U(i) * ( greenMatrix(i,i)*greenMatrix(i+L,i+L) - greenMatrix(i, i+L)*greenMatrix(i+L, i) );
        Venergy += U(i) * (-1.0* cpcpMatrix(i,i+L)*cmcmMatrix(i,i+L) );
    }

    //Add mu and pinning field
    const TensorHao< double, 1> &mu = hubbardSOC->getMu();
    const TensorHao< double, 1> &hx = hubbardSOC->getHx();
    const TensorHao< double, 1> &hy = hubbardSOC->getHy();
    const TensorHao< double, 1> &hz = hubbardSOC->getHz();
    for(size_t i = 0; i < L; ++i)
    {
        Renergy += ( -mu(i) + hz(i)*0.5 ) * greenMatrix(i,i);
        Renergy += ( -mu(i) - hz(i)*0.5 ) * greenMatrix(i+L,i+L);
        Renergy += complex<double>( hx(i)*0.5, -hy(i)*0.5 ) * greenMatrix(i, i+L);
        Renergy += complex<double>( hx(i)*0.5,  hy(i)*0.5 ) * greenMatrix(i+L, i);
    }

    HNum += ( ( Kenergy + Venergy + Renergy ) * denIncrement );
    KNum += ( Kenergy * denIncrement );
    VNum += ( Venergy * denIncrement );
    RNum += ( Renergy * denIncrement );

}
