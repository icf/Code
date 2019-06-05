//
// Created by Hao Shi on 12/18/18.
//

#include "../include/bcs.h"

using namespace std;
using namespace tensor_hao;

Bcs::Bcs() { }

Bcs::~Bcs() { }

void Bcs::run()
{
    initialParameters();

    for(size_t i = 0; i < ( method.annealStep+1 ); ++i)
    {
        findPairing();

        if( variationalEnergy < minimumEnergy )
        {
            minimumEnergy = variationalEnergy;
            minimumStateF = variationalStateF;
            cout<<"The minimumEnergy is updated to "<<fixed<<setprecision(16)<<minimumEnergy<<"\n"<<endl;
            writeFile(minimumEnergy, "minimumEnergy.dat");
            minimumStateF.write("minimumStateF.dat");
            //icf: match to the request of GpBCS input
            //Wf
            TensorHao<complex<double>,2> pBCS_PhiT(2*L,2*L); pBCS_PhiT=complex <double> (0.0,0.0);
            for(size_t temp_i=1-1;temp_i<=L-1;temp_i++){
            for(size_t temp_j=1-1;temp_j<=L-1;temp_j++){
               pBCS_PhiT(temp_i,temp_j+L)=minimumStateF(temp_i,temp_j);
               pBCS_PhiT(temp_i+L,temp_j)=-1.0*pBCS_PhiT(temp_i,temp_j+L);
            }
            }
            pBCS_PhiT.write("pBCS_PhiT.Wf.dat");
            //Orbital
            TensorHao<double,1> tempDia(L);
            BL_NAME(eigen)(minimumStateF,tempDia);
            TensorHao<complex<double>,2> pBCS_PhiT_orbital(2*L,2*L); pBCS_PhiT_orbital=complex <double> (0.0,0.0);
            for(size_t temp_i=1-1;temp_i<=L-1;temp_i++){
            for(size_t temp_j=1-1;temp_j<=L-1;temp_j++){
               pBCS_PhiT_orbital(temp_i,temp_j+L)=minimumStateF(temp_i,temp_j);
               pBCS_PhiT_orbital(temp_i+L,temp_j)=pBCS_PhiT_orbital(temp_i,temp_j+L);
            }
            }
            pBCS_PhiT_orbital.write("pBCS_PhiT.T.dat");
            //Occupancy
            TensorHao<complex<double>,2> pBCS_PhiT_occupancy(2*L,2*L); pBCS_PhiT_occupancy=complex <double> (0.0,0.0);
            for(size_t temp_i=1-1;temp_i<=L-1;temp_i++){
               pBCS_PhiT_occupancy(temp_i,temp_i+L)=tempDia(temp_i);
               pBCS_PhiT_occupancy(temp_i+L,temp_i)=-1.0*pBCS_PhiT_occupancy(temp_i,temp_i+L);
            }
            pBCS_PhiT_occupancy.write("pBCS_PhiT.Gdia.dat");
        }

        if( i!=method.annealStep )
        {
            double anneal = method.annealMagnitude*(method.annealStep-i)/method.annealStep;
            cout<<"Anneal the order parameter: "<<anneal<<"\n"<<endl;
            annealPairing( anneal );
        }
    }

    prepareStop();
}

void Bcs::initialParameters()
{
    if( MPIRank()==0 ) method.read("bcs_param");
    MPIBcast(method);

    randomHaoInit(method.seed, 1);
    if( method.seed != -1 ) randomHaoSave();

    readAndBcastModel("model_param");

    setH0();

    initialHartreeFockEssential();

    minimumEnergy = 1e300;
}

void Bcs::prepareStop()
{
    writeFile(bcsMu, "bcsMu.dat");
    writeFile(bcsN, "bcsN.dat");

    pairing.write("pairing.dat");
}

