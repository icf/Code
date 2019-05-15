//
// Created by boruoshihao on 1/13/17.
// Modefied by icf on 04/29/2019
//

#ifndef AFQMCLAB_HUBBARDSOCBCSSDMEASUREOBSERVE_H
#define AFQMCLAB_HUBBARDSOCBCSSDMEASUREOBSERVE_H

#include "HubbardSOCMeasureCommuteBCSSD.h"

class HubbardSOCMeasureObserveBCSSD : public HubbardSOCMeasureCommuteBCSSD
{
 private:
    tensor_hao::TensorHao<std::complex<double>, 2> greenMatrixNum;
    tensor_hao::TensorHao<std::complex<double>, 2> densityDensityNum;
    tensor_hao::TensorHao<std::complex<double>, 2> splusSminusNum;
    tensor_hao::TensorHao<std::complex<double>, 2> sminusSplusNum;
    tensor_hao::TensorHao<std::complex<double>, 2> spairSpairNum;

 public:
    HubbardSOCMeasureObserveBCSSD();
    HubbardSOCMeasureObserveBCSSD(const HubbardSOC &hubbardSOC_);
    ~HubbardSOCMeasureObserveBCSSD();

    void reSet();
    void addMeasurement(BCSSDOperation &bcssdOperation, std::complex<double> denIncrement);
    void addMeasurementBCSBP(BCSSDOperation &bcssdOperation, std::complex<double> denIncrement);
    const TensorHao<complex<double>, 2> returnGreenMatrix();
    void write() const;
    double getMemory() const;

 private:
    void addGreenMatrix(BCSSDOperation &bcssdOperation, std::complex<double> denIncrement);
    void addGreenMatrixBCSBP(BCSSDOperation &bcssdOperation, std::complex<double> denIncrement);
    //void addDensityDensity(BCSSDOperation &bcssdOperation, std::complex<double> denIncrement);
    //void addSplusSminus(BCSSDOperation &bcssdOperation, std::complex<double> denIncrement);  //icf: we don't do this for now
    //void addSminusSplus(BCSSDOperation &bcssdOperation, std::complex<double> denIncrement);
    //void addSpairSpair(BCSSDOperation &bcssdOperation, std::complex<double> denIncrement);

    HubbardSOCMeasureObserveBCSSD(const HubbardSOCMeasureObserveBCSSD& x);
    HubbardSOCMeasureObserveBCSSD & operator  = (const HubbardSOCMeasureObserveBCSSD& x);
};

#endif //AFQMCLAB_HUBBARDSOCBCSSDMEASUREOBSERVE_H
