//
// Created by Hao Shi on 8/1/17.
//

#include "../include/HubbardSOCMeasureFixedBCSSD.h"
#include "afqmclab.h"

using namespace std;
using namespace tensor_hao;

HubbardSOCMeasureFixedBCSSD::HubbardSOCMeasureFixedBCSSD()
{
    initModelWalkerNullptr();
    reSet();
}

HubbardSOCMeasureFixedBCSSD::HubbardSOCMeasureFixedBCSSD(const HubbardSOC &hubbardSOC_, const BCS &walkerLeft_)
{
    setModelWalker(hubbardSOC_, walkerLeft_);
    reSet();
}

HubbardSOCMeasureFixedBCSSD::~HubbardSOCMeasureFixedBCSSD()
{

}

void HubbardSOCMeasureFixedBCSSD::initModelWalkerNullptr()
{
    hubbardSOC = nullptr;
    walkerLeft = nullptr;
}

void HubbardSOCMeasureFixedBCSSD::setModelWalker(const HubbardSOC &hubbardSOC_, const BCS &walkerLeft_)
{
    if( 2*hubbardSOC_.getL() != walkerLeft_.getL() ) {cout<<"Model L does not consistent with walker L!"<<endl; exit(1);}

    hubbardSOC = &hubbardSOC_;
    walkerLeft = &walkerLeft_;
}
void HubbardSOCMeasureFixedBCSSD::reSet()
{
    complex<double> zero(0,0);
    den = zero;
    KNum = zero; VNum = zero; RNum = zero; HNum = zero;
    NupNum = zero; NdnNum = zero; SplusNum = zero; SminusNum = zero;
    NupTotNum = zero; NdnTotNum = zero; SplusTotNum = zero; SminusTotNum = zero;
}

complex<double> HubbardSOCMeasureFixedBCSSD::returnEnergy()
{
    complex<double> Htot   = MPISum(HNum);
    complex<double> denTot = MPISum(den);
    complex<double> energy;

    if( MPIRank() == 0 ) energy = Htot/denTot;
    MPIBcast(energy);
    return energy;
}

void HubbardSOCMeasureFixedBCSSD::addMeasurement(BCSSDOperation &bcssdOperation, complex<double> denIncrement)
{
    checkWalkerLeft(bcssdOperation);

    den += denIncrement;

    addEnergy(bcssdOperation, denIncrement);
}

NiupNidnForce HubbardSOCMeasureFixedBCSSD::getForce(const NiupNidn &niupNidn, BCSSDOperation &bcssdOperation, double cap)
{
    checkWalkerLeft(bcssdOperation);

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
    else
    {
        cout<<"Error! Can not find the matched decompType! "<<decompType<<endl;
        exit(1);
    }

    NiupNidnForce force(halfL);
    const TensorHao<complex<double>, 1> &gamma = niupNidn.getGamma();
    for (size_t i = 0; i < halfL; ++i)
    {
        force(i) = ( gamma(i) * backGround(i) ).real();
        if( force(i) >  cap ) force(i) =  cap;
        if( force(i) < -cap ) force(i) = -cap;
    }

    return force;
}

void HubbardSOCMeasureFixedBCSSD::write() const
{
    writeThreadSum(den, "den.dat", ios::app);

    writeThreadSum(KNum, "KNum.dat", ios::app);
    writeThreadSum(VNum, "VNum.dat", ios::app);
    writeThreadSum(RNum, "RNum.dat", ios::app);
    writeThreadSum(HNum, "HNum.dat", ios::app);

    writeThreadSum(NupNum.size(),    NupNum.data(),    "NupNum.dat",    ios::app);
    writeThreadSum(NdnNum.size(),    NdnNum.data(),    "NdnNum.dat",    ios::app);
    writeThreadSum(SplusNum.size(),  SplusNum.data(),  "SplusNum.dat",  ios::app);
    writeThreadSum(SminusNum.size(), SminusNum.data(), "SminusNum.dat", ios::app);

    writeThreadSum(NupTotNum,   "NupTotNum.dat",   ios::app);
    writeThreadSum(NdnTotNum,   "NdnTotNum.dat",   ios::app);
    writeThreadSum(SplusTotNum, "SplusTotNum.dat", ios::app);
    writeThreadSum(SminusTotNum, "SminusTotNum.dat", ios::app);
}

double HubbardSOCMeasureFixedBCSSD::getMemory() const
{
    double mem(0.0);
    mem += 8.0;
    mem += 8.0;
    mem += 16.0;
    mem += 16.0*4;
    mem += NupNum.getMemory()+NdnNum.getMemory()+SplusNum.getMemory()+SminusNum.getMemory();
    mem += 16.0*4;
    return mem;
}

HubbardSOCMeasureFixedBCSSD::HubbardSOCMeasureFixedBCSSD(const HubbardSOCMeasureFixedBCSSD &x)
{
}

HubbardSOCMeasureFixedBCSSD & HubbardSOCMeasureFixedBCSSD::operator=(const HubbardSOCMeasureFixedBCSSD &x)
{
    return *this;
}


void HubbardSOCMeasureFixedBCSSD::checkWalkerLeft(const BCSSDOperation &bcssdOperation)
{
    if( walkerLeft != bcssdOperation.getWalkerLeft() )
    {
        cout<<"Error!!! HubbardSOCMeasureFixedBCSSD only accept bcssdOperation with fixed BCS!"<<endl;
        exit(1);
    }
}

void HubbardSOCMeasureFixedBCSSD::addEnergy(BCSSDOperation &bcssdOperation, complex<double> denIncrement)
{
    size_t L = hubbardSOC->getL(); size_t L2 = L*2;  

    const TensorHao< complex<double>, 2 > &greenMatrix = bcssdOperation.returnGreenMatrix();
    const TensorHao<complex<double>, 1> &greenDiagonal = bcssdOperation.returnGreenDiagonal();
    const TensorHao<complex<double>, 1> &greenOffDiagonal = bcssdOperation.returnGreenOffDiagonal();
    const TensorHao< complex<double>, 2 > &cmcmMatrix = bcssdOperation.returnCmCmMatrix();
    const TensorHao< complex<double>, 2 > &cpcpMatrix = bcssdOperation.returnCpCpMatrix();

    complex<double> Kenergy(0,0), Venergy(0,0), Renergy(0,0);

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
        Renergy += ( -mu(i) + hz(i)*0.5 ) * greenDiagonal(i);
        Renergy += ( -mu(i) - hz(i)*0.5 ) * greenDiagonal(i+L);
        Renergy += complex<double>( hx(i)*0.5, -hy(i)*0.5 ) * greenOffDiagonal(i);
        Renergy += complex<double>( hx(i)*0.5,  hy(i)*0.5 ) * greenOffDiagonal(i+L);
    }

    KNum += ( Kenergy * denIncrement );
    VNum += ( Venergy * denIncrement );
    RNum += ( Renergy * denIncrement );
    HNum += ( ( Kenergy + Venergy + Renergy ) * denIncrement );
}


