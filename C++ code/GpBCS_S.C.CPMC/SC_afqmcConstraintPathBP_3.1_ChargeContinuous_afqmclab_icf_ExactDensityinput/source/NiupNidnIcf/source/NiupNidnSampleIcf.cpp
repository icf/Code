//
// Created by boruoshihao on 1/9/17.
//

#include "../include/NiupNidnSampleIcf.h"

using namespace std;
using namespace tensor_hao;

NiupNidnSampleIcf::NiupNidnSampleIcf():logw(0.0) { }

NiupNidnSampleIcf::NiupNidnSampleIcf(size_t L):logw(0.0)
{
    diag00.resize(L);
    diag10.resize(L);
    diag01.resize(L);
    diag11.resize(L);
}

NiupNidnSampleIcf::NiupNidnSampleIcf(const NiupNidnSampleIcf &x) { copy_deep(x);}

NiupNidnSampleIcf::NiupNidnSampleIcf(NiupNidnSampleIcf &&x) { move_deep(x); }

NiupNidnSampleIcf::~NiupNidnSampleIcf() { }

NiupNidnSampleIcf &NiupNidnSampleIcf::operator=(const NiupNidnSampleIcf &x) { copy_deep(x); return *this; }

NiupNidnSampleIcf &NiupNidnSampleIcf::operator=(NiupNidnSampleIcf &&x) { move_deep(x); return *this; }

size_t NiupNidnSampleIcf::getL() const { return diag00.size(); }

double NiupNidnSampleIcf::getMemory() const
{
    return 16.0+diag00.getMemory()+diag01.getMemory()+diag10.getMemory()+diag11.getMemory();
}

void NiupNidnSampleIcf::copy_deep(const NiupNidnSampleIcf &x)
{
    logw = x.logw;
    diag00 = x.diag00;
    diag10 = x.diag10;
    diag01 = x.diag01;
    diag11 = x.diag11;
}

void NiupNidnSampleIcf::move_deep(NiupNidnSampleIcf &x)
{
    logw = x.logw;
    diag00 = move( x.diag00 );
    diag10 = move( x.diag10 );
    diag01 = move( x.diag01 );
    diag11 = move( x.diag11 );
}
