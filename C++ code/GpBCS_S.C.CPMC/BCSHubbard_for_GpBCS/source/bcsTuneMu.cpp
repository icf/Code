//
// Created by Hao Shi on 12/30/18.
//

#include "../include/bcs.h"

using namespace std;
using namespace tensor_hao;

void Bcs::findBcsMu()
{
    findBcsMuPlusBcsMuMinus();
    findBcsMuBetweenPlusAndMinus();
}

void Bcs::calculateBcsNFromPairingAndMu()
{
    setMeanFieldVectorsValuesFromPairingAndMu();
    setVariationalStateFromMeanFieldVectors();
    calculateVariationalStateFGreen();
    calculateNparticle();
}

void Bcs::findBcsMuPlusBcsMuMinus()
{
    calculateBcsNFromPairingAndMu();

    cout<<setprecision(8)<<endl;

    if( bcsN < N )
    {
        while ( bcsN<N )
        {
            bcsMu += method.deltaMu;
            cout<<"  Looking for bcsMuPlus, N particle is "<<bcsN<<", bcsMu "<<bcsMu<<"."<<endl;
            calculateBcsNFromPairingAndMu();
        }

        bcsMuPlus  = bcsMu + 0.001; //avoid bcsN==N
        bcsMuMinus = bcsMu - method.deltaMu;

    }
    else if ( bcsN>N )
    {
        while ( bcsN>N )
        {
            bcsMu -= method.deltaMu;
            cout<<"  Looking for bcsMuMinus, N particle is "<<bcsN<<", bcsMu "<<bcsMu<<"."<<endl;
            calculateBcsNFromPairingAndMu();
        }
        bcsMuPlus  = bcsMu + method.deltaMu;
        bcsMuMinus = bcsMu - 0.001; //avoid bcsN==N
    }
    else
    {
        bcsMuPlus  = bcsMu + method.deltaMu;
        bcsMuMinus = bcsMu - method.deltaMu;
    }
    cout<<"  Find bcsMuMinus and bcsMuPlus: "<<bcsMuMinus<<" "<<bcsMuPlus<<"\n"<<endl;
}

void Bcs::findBcsMuBetweenPlusAndMinus()
{
    cout<<setprecision(8)<<endl;

    while ( abs(bcsN - N) > method.Neta )
    {
        cout<<"  Looking for bcsMu, N particle is "<<bcsN<<", bcsMu "<<bcsMu<<"."<<endl;

        bcsMu = ( bcsMuPlus + bcsMuMinus )/2.0;

        calculateBcsNFromPairingAndMu();

        if( bcsN < N )
        {
            bcsMuMinus = bcsMu;
        }
        else
        {
            bcsMuPlus = bcsMu;
        }
    }

    cout<<"  Find bcsMu and bcsN: "<<bcsMu<<" "<<bcsN<<" "<<N<<"\n\n\n"<<endl;
}