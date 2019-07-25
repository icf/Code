//
// Created by boruoshihao on 1/8/17.
//

#ifndef AFQMCLAB_NIUPNIDNICF_H
#define AFQMCLAB_NIUPNIDNICF_H

#include "NiupNidnAuxIcf.h"
#include "NiupNidnForceIcf.h"
#include "NiupNidnSampleIcf.h"

class NiupNidn
{
 private:
    size_t L;
    std::string decompType;  //densityCharge, densitySpin, hopCharge, hopSpin, pairCharge, pairSpin
    double dtUSum;
    tensor_hao::TensorHao<double,1> dtU;
    tensor_hao::TensorHao<std::complex<double>,1> gamma;
    tensor_hao::TensorHao<std::complex<double>,1> constDiag00, constDiag10, constDiag01, constDiag11;

 public:
    NiupNidn();
    NiupNidn(double dt,
             const std::string &decompTypeIn,
             const tensor_hao::TensorHao<double, 1> &U,
             const tensor_hao::TensorHao<double, 1> &mu,
             const tensor_hao::TensorHao<double, 1> &hx,
             const tensor_hao::TensorHao<double, 1> &hy,
             const tensor_hao::TensorHao<double, 1> &hz);
    NiupNidn(const NiupNidn& x);
    NiupNidn(NiupNidn&& x);
    ~NiupNidn();

    NiupNidn & operator  = (const NiupNidn& x);
    NiupNidn & operator  = (NiupNidn&& x);

    size_t getL() const;
    const std::string &getDecompType() const;
    double getDtUSum() const;
    const tensor_hao::TensorHao<double, 1> &getDtU() const;
    const tensor_hao::TensorHao<std::complex<double>, 1> &getGamma() const;
    const tensor_hao::TensorHao<std::complex<double>, 1> &getConstDiag00() const;
    const tensor_hao::TensorHao<std::complex<double>, 1> &getConstDiag10() const;
    const tensor_hao::TensorHao<std::complex<double>, 1> &getConstDiag01() const;
    const tensor_hao::TensorHao<std::complex<double>, 1> &getConstDiag11() const;

    NiupNidnForceIcf readForce(const std::string &filename) const;
    NiupNidnAuxIcf sampleAuxFromForce(const NiupNidnForceIcf &force) const;
    double logProbOfAuxFromForce(const NiupNidnAuxIcf &aux, const NiupNidnForceIcf &force) const;
    NiupNidnSampleIcf getTwoBodySampleFromAux(const NiupNidnAuxIcf &aux) const;
    NiupNidnSampleIcf getTwoBodySampleFromAuxForce(const NiupNidnAuxIcf &aux, const NiupNidnForceIcf &force) const;

    size_t getAuxSize() const;
    size_t getAuxDiffSize(const NiupNidnAuxIcf &auxOne, const NiupNidnAuxIcf &auxTwo) const;
    double getMemory() const;

 private:
    void copy_deep(const NiupNidn &x);
    void move_deep(NiupNidn &x);

    void setGamma();
    void setConstDiag(double dt,
                      const tensor_hao::TensorHao<double,1> &U,
                      const tensor_hao::TensorHao<double,1> &mu,
                      const tensor_hao::TensorHao<double,1> &hx,
                      const tensor_hao::TensorHao<double,1> &hy,
                      const tensor_hao::TensorHao<double,1> &hz);
    void setTwoBodySampleMatrix(NiupNidnSampleIcf &twoBodySample, const NiupNidnAuxIcf &aux) const;
};

#endif //AFQMCLAB_NIUPNIDNICF_H
