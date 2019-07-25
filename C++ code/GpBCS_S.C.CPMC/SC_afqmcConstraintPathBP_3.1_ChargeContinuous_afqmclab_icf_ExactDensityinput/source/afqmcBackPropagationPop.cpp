//
// Created by Hao Shi on 2/4/19.
//

#include "../include/afqmcBackPropagationPop.h"

using namespace std;

AfqmcBackPropagationPop::AfqmcBackPropagationPop()
{

}

AfqmcBackPropagationPop::AfqmcBackPropagationPop(size_t BPNumber_in,
                                                 vector<TwoBodyAux> &twoBodyAuxBPMeasure_in,
                                                 vector<complex<double>> &reWeightBPMeasure_in,
                                                 WalkerRight &walkerBPMeasure_in)
{
    BPNumber = BPNumber_in;
    twoBodyAuxBPMeasure = &twoBodyAuxBPMeasure_in;
    reWeightBPMeasure = &reWeightBPMeasure_in;
    walkerBPMeasure = &walkerBPMeasure_in;
    setNbuf();
}

AfqmcBackPropagationPop::~AfqmcBackPropagationPop()
{

}

int AfqmcBackPropagationPop::getNbuf() const
{
    return Nbuf;
}

AfqmcBackPropagationPop &AfqmcBackPropagationPop::operator=(const AfqmcBackPropagationPop &x)
{
    for (size_t i = 0; i < BPNumber; ++i)
    {
        twoBodyAuxBPMeasure->operator[](i) = x.twoBodyAuxBPMeasure->operator[](i);
        reWeightBPMeasure->operator[](i) = x.reWeightBPMeasure->operator[](i);
    }
    *walkerBPMeasure = *(x.walkerBPMeasure);
    return *this;
}

#ifdef MPI_HAO

vector<char> AfqmcBackPropagationPop::pack() const
{
    vector<char> buf(Nbuf);

    int posit=0;
    for (size_t i = 0; i < BPNumber; ++i)
    {
        TwoBodyAux &twoBodyAux = twoBodyAuxBPMeasure->operator[](i);
        MPI_Pack( twoBodyAux.data(), twoBodyAux.size(), MPI_DOUBLE_COMPLEX, buf.data(), buf.size(), &posit, MPI_COMM_WORLD);
        complex<double> &reWeight = reWeightBPMeasure->operator[](i);
        MPI_Pack(&reWeight, 1, MPI_DOUBLE_COMPLEX, buf.data(), buf.size(), &posit, MPI_COMM_WORLD);
    }
    walkerBPMeasure->pack( buf, posit );

    if(posit!=Nbuf) {cout<<"ERROR in pack!!! posit does not equal Nbuf! "<<posit<<" "<<Nbuf<<endl; exit(1);}

    return buf;
}

void AfqmcBackPropagationPop::unpack(const vector<char>& buf)
{
    int bufSize=buf.size();
    if(bufSize!=Nbuf) {cout<<"ERROR in unpack!!! Size of input buf does not equal Nbuf! "<<endl; exit(1);}

    int posit=0;
    for (size_t i = 0; i < BPNumber; ++i)
    {
        TwoBodyAux &twoBodyAux = twoBodyAuxBPMeasure->operator[](i);
        MPI_Unpack(buf.data(), buf.size(), &posit, twoBodyAux.data(), twoBodyAux.size(), MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
        complex<double> &reWeight = reWeightBPMeasure->operator[](i);
        MPI_Unpack(buf.data(), buf.size(), &posit, &reWeight, 1, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
    }
    walkerBPMeasure->unpack( buf, posit );

    if(posit!=Nbuf) {cout<<"ERROR in unpack!!! posit does not equal Nbuf! "<<posit<<" "<<Nbuf<<endl; exit(1);}
}

#endif

void AfqmcBackPropagationPop::setNbuf()
{
    Nbuf  = 0;
    Nbuf += BPNumber * sizeof( complex<double> ) * twoBodyAuxBPMeasure->operator[](0).size();
    Nbuf += BPNumber * sizeof( complex<double> );
    Nbuf += walkerBPMeasure->returnNbuf();
}
