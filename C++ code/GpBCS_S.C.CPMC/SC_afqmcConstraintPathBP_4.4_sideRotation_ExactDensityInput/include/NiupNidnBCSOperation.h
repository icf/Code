//
// Created by boruoshihao on 1/10/17.
// Modefied by icf on 05/06/2019
//

#ifndef AFQMCLAB_NIUPNIDNBCSOPERATION_H
#define AFQMCLAB_NIUPNIDNBCSOPERATION_H

#include <tuple>
#include "BCS.h"
#include "afqmclab.h"

class NiupNidnSampleBCSOperation
{
 public:
    NiupNidnSampleBCSOperation();
    ~NiupNidnSampleBCSOperation();

    void applyToRight(const NiupNidnSample &oneBody, const BCS &walker, BCS &walkerNew) const;
    void applyToLeft(const NiupNidnSample &oneBody, const BCS &walker, BCS &walkerNew) const;
};

#endif //AFQMCLAB_NIUPNIDNBCSOPERATION_H
