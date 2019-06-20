//
// Created by Hao Shi on 8/1/17.
// Modefied by icf on 04/29/2019
//

#ifndef AFQMCLAB_HUBBARDSOCMEASUREFIXEDBCSSD_H
#define AFQMCLAB_HUBBARDSOCMEASUREFIXEDBCSSD_H

#include "afqmclab.h"   
//#include "HubbardSOC.h"
#include "BCSSDOperation.h"

class HubbardSOCMeasureFixedBCSSD
{
 private:
    const HubbardSOC * hubbardSOC;
    const BCS *walkerLeft;

    std::complex<double> den;
    std::complex<double> KNum, VNum, RNum, HNum;
    tensor_hao::TensorHao<std::complex<double>, 1> NupNum, NdnNum, SplusNum, SminusNum;
    std::complex<double> NupTotNum, NdnTotNum, SplusTotNum, SminusTotNum;

    //tensor_hao::TensorHao<std::complex<double>,2> wfDaggerK;   //icf: no need of this

 public:
    HubbardSOCMeasureFixedBCSSD();
    HubbardSOCMeasureFixedBCSSD(const HubbardSOC& hubbardSOC_, const BCS& walkerLeft_);
    ~HubbardSOCMeasureFixedBCSSD();

    void initModelWalkerNullptr();
    void setModelWalker(const HubbardSOC& hubbardSOC_, const BCS& walkerLeft_);
    void reSet();
    std::complex<double> returnEnergy();
    void addMeasurement(BCSSDOperation &bcssdOperation, std::complex<double> denIncrement);
    NiupNidnForce getForce(const NiupNidn &niupNidn, BCSSDOperation &bcssdOperation, double cap=1.0);

    void write() const;
    double getMemory() const;

 private:
    HubbardSOCMeasureFixedBCSSD(const HubbardSOCMeasureFixedBCSSD& x);
    HubbardSOCMeasureFixedBCSSD & operator = (const HubbardSOCMeasureFixedBCSSD& x);

    //void initWfDaggerK();
    void checkWalkerLeft(const BCSSDOperation &bcssdOperation);
    void addEnergy(BCSSDOperation &bcssdOperation, std::complex<double> denIncrement);
    //void addNupNdn(BCSSDOperation &bcssdOperation, std::complex<double> denIncrement);
    //void addSplusSminus(BCSSDOperation &bcssdOperation, std::complex<double> denIncrement); //icf: we do not do it for now
};


#endif //AFQMCLAB_HUBBARDSOCMEASUREFIXEDBCSSD_H
