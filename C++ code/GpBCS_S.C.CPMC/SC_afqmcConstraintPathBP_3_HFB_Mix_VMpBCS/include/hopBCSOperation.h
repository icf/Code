//
// Created by boruoshihao on 12/28/16.
// Modefied by icf on 05/06/2019
//

#ifndef AFQMCLAB_HOPBCSOPERATION_H
#define AFQMCLAB_HOPBCSOPERATION_H

#include "BCS.h"
#include "afqmclab.h"

class HopBCSOperation
{
 public:
    HopBCSOperation();
    ~HopBCSOperation();

    void applyToRight(const Hop &oneBody, const BCS &walker, BCS &walkerNew) const;
    void applyToLeft(const Hop &oneBody, const BCS &walker, BCS &walkerNew) const;

 private:
    void checkAndResize(const Hop &oneBody, const BCS &walker, BCS &walkerNew) const;
};

#endif //AFQMCLAB_HOPBCSOPERATION_H
