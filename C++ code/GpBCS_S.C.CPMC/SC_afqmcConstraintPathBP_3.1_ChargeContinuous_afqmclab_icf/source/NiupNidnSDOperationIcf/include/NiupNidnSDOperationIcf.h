//
// Created by boruoshihao on 1/10/17.
//

#ifndef AFQMCLAB_NIUPNIDNSDOPERATIONICF_H
#define AFQMCLAB_NIUPNIDNSDOPERATIONICF_H

#include <tuple>
//#include "../../../walker/SD/include/SD.h"
//#include "../../../walkerWalkerOperation/SDSDOperation/include/SDSDOperation.h"
#include "afqmclab.h"
#include "NiupNidnIcf.h"

class NiupNidnSampleSDOperationIcf
{
 public:
    NiupNidnSampleSDOperationIcf();
    ~NiupNidnSampleSDOperationIcf();

    void applyToRight(const NiupNidnSample &oneBody, const SD &walker, SD &walkerNew) const;
    void applyToLeft(const NiupNidnSample &oneBody, const SD &walker, SD &walkerNew) const;
};

#endif //AFQMCLAB_NIUPNIDNSDOPERATIONIcf_H
