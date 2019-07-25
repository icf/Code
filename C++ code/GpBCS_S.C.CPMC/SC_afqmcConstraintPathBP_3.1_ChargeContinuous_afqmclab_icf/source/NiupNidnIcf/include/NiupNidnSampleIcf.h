//
// Created by boruoshihao on 1/9/17.
//

#ifndef AFQMCLAB_NIUPNIDNSAMPLEICF_H
#define AFQMCLAB_NIUPNIDNSAMPLEICF_H

#include "afqmclab.h"

class NiupNidnSampleIcf
{
 public:
    std::complex<double> logw;
    tensor_hao::TensorHao<std::complex<double>,1> diag00, diag10, diag01, diag11;

    NiupNidnSampleIcf();
    NiupNidnSampleIcf(size_t L);
    NiupNidnSampleIcf(const NiupNidnSampleIcf &x);
    NiupNidnSampleIcf(NiupNidnSampleIcf &&x);
    ~NiupNidnSampleIcf();

    NiupNidnSampleIcf & operator  = (const NiupNidnSampleIcf& x);
    NiupNidnSampleIcf & operator  = (NiupNidnSampleIcf&& x);

    size_t getL() const;
    double getMemory() const;

 private:
    void copy_deep(const NiupNidnSampleIcf &x);
    void move_deep(NiupNidnSampleIcf &x);
};

#endif //AFQMCLAB_NIUPNIDNSAMPLEICF_H
