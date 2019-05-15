//
// Created by boruoshihao on 1/13/17.
// Modefied by icf on 04/29/2019
//

#ifndef AFQMCLAB_HUBBARDSOCBCSSDMEASURECOMMUTE_H
#define AFQMCLAB_HUBBARDSOCBCSSDMEASURECOMMUTE_H

#include "afqmclab.h"
#include "BCSSDOperation.h"

class HubbardSOCMeasureCommuteBCSSD
{
 private:
    const HubbardSOC * hubbardSOC;
    std::complex<double> HNum, KNum, VNum, RNum;

 public:

    std::complex<double> den;

    HubbardSOCMeasureCommuteBCSSD();
    HubbardSOCMeasureCommuteBCSSD(const HubbardSOC &hubbardSOC_);
    ~HubbardSOCMeasureCommuteBCSSD();

    const HubbardSOC *getHubbardSOC() const;

    void initModelNullptr();
    void setModel(const HubbardSOC &hubbardSOC_);
    void reSet();
    std::complex<double> returnEnergy();
    void addMeasurement(BCSSDOperation &bcssdOperation, std::complex<double> denIncrement);
    void addMeasurementBCSBP(BCSSDOperation &bcssdOperation, std::complex<double> denIncrement);
    NiupNidnForce getForce(const NiupNidn &niupNidn, BCSSDOperation &bcssdOperation, double cap=1.0);

    void write() const;
    void writeKNumVumRum() const;
    double getMemory() const;

 private:
    HubbardSOCMeasureCommuteBCSSD(const HubbardSOCMeasureCommuteBCSSD& x);
    HubbardSOCMeasureCommuteBCSSD & operator  = (const HubbardSOCMeasureCommuteBCSSD& x);
    void addEnergy(const tensor_hao::TensorHao<std::complex <double>, 2> &greenMatrix, const tensor_hao::TensorHao<std::complex <double>, 2> &cmcmMatrix, const tensor_hao::TensorHao<std::complex <double>, 2> &cpcpMatrix, std::complex <double> denIncrement);
};

#endif //AFQMCLAB_HUBBARDSOCBCSSDMEASURECOMMUTE_H
