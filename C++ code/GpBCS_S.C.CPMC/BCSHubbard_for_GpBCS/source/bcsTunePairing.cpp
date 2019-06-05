//
// Created by Hao Shi on 12/30/18.
//

#include "../include/bcs.h"

using namespace std;
using namespace tensor_hao;

void Bcs::findPairing()
{
    size_t i;
    for(i = 0; i < method.maxIterateStep; ++i)
    {
        findBcsMu();
        calculatePairing();
        calculateVariationalEnergy();
        cout<<"VariationalEnergy is "<<variationalEnergy<<endl;
        cout<<pairing<<endl;

        calculateConvergeValue();
        if( convergeValue < method.convergeTolerance ) break;

        relaxPairing();
    }
    if( i < method.maxIterateStep ) cout<<"\nIt take "<<i<<" step to converge.\n"<<endl;
    else cout<<"\nWARNING!!! Does not converge!!! \n"<<endl;
    cout<<"ConvergeValue is  "<<convergeValue<<" \n"<<endl;
}

void Bcs::calculateConvergeValue()
{
//    double difference(0), total(0);
//
//    for(size_t i = 0; i < L; ++i)
//    {
//        difference += abs( pairing(i) - pairingBefore(i) );
//        total += abs( pairing(i) );
//    }
//    convergeValue = ( difference * L )/ total;

    convergeValue = abs( (variationalEnergy-variationalEnergyBefore)/variationalEnergy );
}

void Bcs::relaxPairing()
{
    for(size_t i = 0; i < L; ++i)
    {
        pairing(i) = pairingBefore(i) + method.relaxMagnitude * ( pairing(i)-pairingBefore(i) );
    }
}

void Bcs::annealPairing(double anneal)
{
    for(size_t i = 0; i < L; ++i)
    {
        pairing(i) *= (1.0 + (2.0 * uniformHao() - 1.0) * anneal);
    }
}

void Bcs::checkPairing()
{
    double eta = 1e-5;

    for(size_t i = 0; i < L; ++i)
    {
        if( abs( pairing(i).imag() ) > eta )
        {
            cout<<"Error!!! pairing(i) has imaginary part: " << i << " " << pairing(i) <<endl;
            exit(1);
        }
    }
}