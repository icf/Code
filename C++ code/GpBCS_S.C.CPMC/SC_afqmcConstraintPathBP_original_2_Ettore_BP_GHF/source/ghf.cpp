//
// Created by boruoshihao on 4/30/17.
//

#include "../include/ghf.h"

using namespace std;
using namespace tensor_hao;

Ghf::Ghf() { }

Ghf::~Ghf() { }

TensorHao< complex<double>, 2 > Ghf::run()
{

    initialParameters();

    for(size_t i = 0; i < ( method.annealStep+1 ); ++i)
    {
        selfConsistentLoop();

        if( variationalEnergy < minimumEnergy )
        {
            minimumEnergy = variationalEnergy;
            minimumState = variationalState;
            cout<<"The minimumEnergy is updated to "<<fixed<<setprecision(16)<<minimumEnergy<<"\n"<<endl;
            writeFile(minimumEnergy, "minimumEnergy.dat");
            minimumState.write("minimumState.dat");
        }

        if( i!=method.annealStep )
        {
            double anneal = method.annealMagnitude*(method.annealStep-i)/method.annealStep;
            cout<<"Anneal the order parameter: "<<anneal<<"\n"<<endl;
            annealOrderParameter( anneal );
        }
    }

    prepareStop();

    return meanFieldVectors;
}

void Ghf::initialParameters()
{
    method.read("ghf_param");

    randomHaoInit(method.seed, 1);
    if( method.seed != -1 ) randomHaoSave();

    model.read("model_param");

    setH0();

    initialHartreeFockEssential();

    minimumEnergy = 1e300;
}

void Ghf::selfConsistentLoop()
{
    size_t i;
    for(i = 0; i < method.maxIterateStep; ++i)
    {
        oneSelfConsistentStep();
        if( convergeValue < method.convergeTolerance ) break;
    }
    if( i < method.maxIterateStep ) cout<<"\nIt take "<<i<<" step to converge.\n"<<endl;
    else cout<<"\nWARNING!!! Does not converge!!! \n"<<endl;
}

void Ghf::oneSelfConsistentStep()
{
    double variationalEnergyBefore= variationalEnergy;
    setMeanFieldVectorsValuesFromOrderParameter();
    setVariationalStateFromMeanFieldVectors();
    backupAndUpdateOrderParameterFromVariationalState();

    setVariationalEnergyFromOrderParameterAndMeanFieldValues();
    relaxOrderParameter();
    if( method.convergeType == "energy" )
    {
        double difference = abs( variationalEnergy-variationalEnergyBefore );
        convergeValue = difference / abs( variationalEnergy );
    }
    else if( method.convergeType == "orderParameter" )
    {
        size_t L  = model.getL();

        double difference(0), total(0);
        for(size_t i = 0; i < 2*L; ++i)
        {
            difference += abs( density(i) - densityBefore(i) );
            difference += abs( spinOrbit(i) - spinOrbitBefore(i) );
            total += abs( density(i) );
            total += abs( spinOrbit(i) );
        }
        convergeValue = difference / total;
    }
}

void Ghf::prepareStop()
{
    SDSDOperation sdsdOperation(minimumState, minimumState);
    HubbardSOCMeasureObserveSDSD meas(model);
    meas.addMeasurement(sdsdOperation, 1.0);
    meas.write();
    //Hao style
    meanFieldVectors.write("Hao_Full_GHF.inputdat");
}
