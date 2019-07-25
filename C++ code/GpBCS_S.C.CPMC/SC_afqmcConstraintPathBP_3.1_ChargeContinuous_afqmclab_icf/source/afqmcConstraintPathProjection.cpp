//
// Created by boruoshihao on 4/16/17.
// Modefied by icf 6/13/2019
//

#include "../include/afqmcConstraintPath.h"
#include "../include/afqmcWalkerPop.h"

using namespace std;
using namespace tensor_hao;

void AfqmcConstraintPath::projectExpHalfDtK()
{
    WalkerRight walkerTemp;
    for(int i = 0; i < method.walkerSizePerThread; ++i)
    {
        if( walkerIsAlive[i] )
        {
            oneBodyWalkerRightOperation.applyToRight(expHalfDtK, walker[i], walkerTemp);
            walker[i] = move( walkerTemp );
        }
    }
}

void AfqmcConstraintPath::projectExpMinusHalfDtK()
{
    WalkerRight walkerTemp;
    for(int i = 0; i < method.walkerSizePerThread; ++i)
    {
        if( walkerIsAlive[i] )
        {
            oneBodyWalkerRightOperation.applyToRight(expMinusHalfDtK, walker[i], walkerTemp); //icf: e^O|i> the def of oneBodyWalkerRightOperation should be modefied if we use different e^O and |i>
            walker[i] = move( walkerTemp );
        }
    }
}

void AfqmcConstraintPath::projectExpMinusDtKExpMinusDtV()
{
    complex<double> logOverlap;
    double logDiff;
    double phase1,phase2;
    WalkerRight walkerTemp;
    complex<double> auxForce;
    for(int i = 0; i < method.walkerSizePerThread; ++i)
    {
        if( walkerIsAlive[i] )
        {
            walkerWalkerOperation.set(phiT, walker[i]);  // get the ovelap <phiT|phi>  

            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//icf: phaseless
            logOverlap = walkerWalkerOperation.returnLogOverlap();
            //Cap weight real
            logDiff = ( logOverlap - logOverlapBackup[i] ).real();
            if( (logDiff-method.logEnergyCap)>1e-12  )
            {
                cout<<"Cap weight in projection: "
                    <<setw(25)<<i+MPIRank()*method.walkerSizePerThread
                    <<setw(25)<<logDiff<<setw(25)<<method.logEnergyCap<<endl;
                walker[i].addLogw( method.logEnergyCap  - logDiff );
                logOverlap += ( method.logEnergyCap - logDiff );
            }

            //Force weight to be real
            phase1 = logOverlap.imag();   
            if (cos(phase1) <= 0.0) { walkerIsAlive[i] = false; continue; }
            walker[i].addLogw( log(cos(phase1)) - complex<double>(0, phase1) );    //icf: Phase cancelling

            logOverlapBackup[i] = logOverlap + log(cos(phase1)) - complex<double>(0, phase1);  //icf: What does logOverlapBackup mean?
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
            if( method.forceType == "constForce" )
            {
                twoBodyAux = expMinusDtV.sampleAuxFromForce(constForce);
                twoBodySample = expMinusDtV.getTwoBodySampleFromAuxForce(twoBodyAux, constForce);
                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//icf: phaseless
                //auxForce = 0.0; for(size_t k = 0; k < model.getCholeskyNumber() ; ++k) auxForce += twoBodyAux(k)*constForce(k);   
                auxForce = 0.0; for(size_t k = 0; k < model.getL() ; ++k) auxForce += twoBodyAux(k)*constForce(k);  
                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
            }
            else if( method.forceType == "dynamicForce" )
            {
                dynamicForce = observeMeasure.getForce(expMinusDtV, walkerWalkerOperation, method.forceCap);
                twoBodyAux = expMinusDtV.sampleAuxFromForce(dynamicForce);
                twoBodySample = expMinusDtV.getTwoBodySampleFromAuxForce(twoBodyAux, dynamicForce);
                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//icf: phaseless
                //auxForce = 0.0; for(size_t k = 0; k < model.getCholeskyNumber() ; ++k) auxForce += twoBodyAux(k)*dynamicForce(k);  //icf: aux*force  for Hubbard "getCholeskyNumber"=Nsite?
                auxForce = 0.0; for(size_t k = 0; k < model.getL() ; ++k) auxForce += -1.0*twoBodyAux(k)*dynamicForce(k);//+dynamicForce(k)*dynamicForce(k);   //Attention: should be zero for half-filling?
                //auxForce=twoBodySample.logw;
                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
            }
            else
            {
                cout<<"Error!!! Do not know method.forceType "<<method.forceType<<endl;
                exit(1);
            }


            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//icf: phaseless
            //Cos projection
            phase2 = auxForce.imag();
            if (cos(phase2) <= 0.0) { walkerIsAlive[i] = false; continue; }
            walker[i].addLogw( log( cos(phase2) ) );
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//


            twoBodySampleWalkerRightOperation.applyToRight(twoBodySample, walker[i], walkerTemp);
            walkerTemp.addLogw( method.dt*method.ET );
            oneBodyWalkerRightOperation.applyToRight(expMinusDtK, walkerTemp, walker[i]);   
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//icf: phaseless
            if( isBP )
            {
                twoBodyAuxBackup[i][BPIndex]=move(twoBodyAux);
                reWeightBackup[i][BPIndex] = log(cos(phase1)) - complex<double>(0, phase1) + log( cos(phase2) );
            }
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
        }
    }

    if( isBP ) tablePerThreadBackup[BPIndex].resize(0);
}

void AfqmcConstraintPath::projectOneStep(size_t &mgsIndex, size_t &popControlIndex)
{
    projectExpMinusDtKExpMinusDtV(); //  e^-0.5T|phiT> --> e^-T e^-V e^-0.5T|phiT>; e^-0.5T from the forward projection before. 

    mgsIndex++;
    if (mgsIndex == method.mgsStep)
    {
        modifyGM();
        mgsIndex = 0;
    }

    popControlIndex++;
    if (popControlIndex == method.popControlStep)
    {
        popControl();
        popControlIndex = 0;
    }

    if( isBP ) BPIndex++;
}

void AfqmcConstraintPath::modifyGM()
{
    for(int i = 0; i < method.walkerSizePerThread; ++i)
    {
        if( walkerIsAlive[i] ) walker[i].stabilize();
    }
}

void AfqmcConstraintPath::popControl()
{

    vector<double> weightPerThread( method.walkerSizePerThread );
    complex<double> logOverlap; double overlapReal;
    double logDiff;
    for(int i = 0; i < method.walkerSizePerThread; ++i)
    {
        if( walkerIsAlive[i] )
        {
            walkerWalkerOperation.set(phiT, walker[i]);

            logOverlap = walkerWalkerOperation.returnLogOverlap();
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
            //Cap weight real
            if( logOverlapBackup.size() > 0 )
            {
                logDiff = ( logOverlap - logOverlapBackup[i] ).real();
                if( (logDiff-method.logEnergyCap)>1e-12  )
                {
                    cout<<"Cap weight in projection before popControl: "
                        <<setw(25)<<i+MPIRank()*method.walkerSizePerThread
                        <<setw(25)<<logDiff<<setw(25)<<method.logEnergyCap<<endl;
                    walker[i].addLogw( method.logEnergyCap  - logDiff );
                    logOverlap += ( method.logEnergyCap - logDiff );
                }
            }

            overlapReal = (exp(logOverlap)).real();
            if (overlapReal <= 0.0) { weightPerThread[i] = 0.0; walkerIsAlive[i] = false; }
            else { weightPerThread[i] = overlapReal; walker[i].addLogw( -logOverlap );}

        }
        else
        {
            weightPerThread[i] = 0.0;
        }
    }

    checkAndResetWalkerIsAlive();
    resetLogOverlapBackup();

    vector<double> weight;

#ifdef MPI_HAO
    if( MPIRank()==0 ) weight.resize( method.walkerSize );
    MPI_Gather( weightPerThread.data(), method.walkerSizePerThread, MPI_DOUBLE_PRECISION,
                weight.data(), method.walkerSizePerThread, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD );
#else
    weight =  move(weightPerThread);
#endif

    if( MPIRank() == 0 ) popCheck(weight);

    vector<int> table;

    if( MPIRank()==0 ) table=popConfiguration( MPISize(), weight );

    vector<AfqmcWalkerPop> walkerPop;
    walkerPop.reserve( method.walkerSizePerThread );
    for(int i=0; i<method.walkerSizePerThread; i++) walkerPop.push_back( AfqmcWalkerPop(walker[i]) );
    populationControl(walkerPop, table);

    if( MPIRank()==0 ) cout<<endl;

    if( isBP )
    {
#ifdef MPI_HAO
        tablePerThreadBackup[BPIndex].resize(method.walkerSizePerThread);
        MPI_Scatter(table.data(), method.walkerSizePerThread, MPI_INT, tablePerThreadBackup[BPIndex].data(), method.walkerSizePerThread, MPI_INT, 0, MPI_COMM_WORLD);
#else
        tablePerThreadBackup[BPIndex] = move(table);
#endif
    }

}

void AfqmcConstraintPath::checkAndResetWalkerIsAlive()
{
    size_t aliveWalkerPerThread(0);
    for (int i = 0; i < method.walkerSizePerThread; ++i)
    {
        if( walkerIsAlive[i] ) aliveWalkerPerThread++;
    }

    size_t aliveWalker = MPISum(aliveWalkerPerThread);

    if( MPIRank()==0 )
    {
        cout<<"Total number of walker is "<<method.walkerSize<<"."<<endl;
        cout<<"Currently "<<aliveWalker<<" walkers are still alive."<<endl;
        cout<<"Currently "<<method.walkerSize-aliveWalker<<" walkers are killed."<<endl;
    }

    for (int i = 0; i < method.walkerSizePerThread ; ++i)
    {
        walkerIsAlive[i] = true;
    }

}

void AfqmcConstraintPath::resetLogOverlapBackup()
{
    logOverlapBackup.resize(method.walkerSizePerThread);
    for (int i = 0; i < method.walkerSizePerThread ; ++i)
    {
        logOverlapBackup[i] = 0.0;
    }
}
