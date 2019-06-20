//
// Created by boruoshihao on 4/16/17.
// Modefied by icf on 05/02/2019
//

#ifndef AFQMCLAB_AFQMCCONSTRAINTPATHDEFINE_H
#define AFQMCLAB_AFQMCCONSTRAINTPATHDEFINE_H

#include "afqmclab.h"
#include "BCSXiao.h"

typedef Hop OneBody;

typedef NiupNidn       TwoBody;
typedef NiupNidnAux    TwoBodyAux;
typedef NiupNidnForce  TwoBodyForce;
typedef NiupNidnSample TwoBodySample;

typedef HubbardSOC Model;

typedef BCS  WalkerLeft;     //change it to BCS
typedef SD  WalkerRight;
typedef HopSDOperation OneBodyWalkerRightOperation;     
typedef HopBCSOperation OneBodyWalkerLeftOperation;       //change it to BCS
typedef NiupNidnSampleSDOperation TwoBodySampleWalkerRightOperation;
typedef NiupNidnSampleBCSOperation TwoBodySampleWalkerLeftOperation;      //change it to BCS
typedef BCSSDOperation WalkerWalkerOperation;     //change it to BCS
typedef HubbardSOCMeasureObserveBCSSD ModelMeasureObserve;      //change it to BCS

#endif //AFQMCLAB_AFQMCCONSTRAINTPATHDEFINE_H
