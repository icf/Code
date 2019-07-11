//
// Created by boruoshihao on 4/16/17.
// Modefied by icf on 05/02/2019
//
#include "../include/afqmcConstraintPath.h"

using namespace std;
using namespace tensor_hao;

void AfqmcConstraintPath::addMixedMeasurement()
{
    complex<double> logOverlap;
    double logDiff;
    WalkerRight walkerTemp;

    for(int i = 0; i < method.walkerSizePerThread; ++i)
    {
        if( walkerIsAlive[i] )
        {
            oneBodyWalkerRightOperation.applyToRight(expHalfDtK, walker[i], walkerTemp);

            walkerWalkerOperation.set(phiT, walkerTemp);

            logOverlap = walkerWalkerOperation.returnLogOverlap();

            //Cap weight real part
            logDiff = ( logOverlap - logOverlapBackup[i] ).real();
            if( (logDiff-method.logEnergyCap)>1e-12 )
            {
                cout<<"Cap weight in measurement: "
                    <<setw(25)<<i+MPIRank()*method.walkerSizePerThread
                    <<setw(25)<<logDiff<<setw(25)<<method.logEnergyCap<<endl;
                logOverlap+= ( method.logEnergyCap - logDiff );
            }

            observeMeasure.addMeasurement(walkerWalkerOperation, exp(logOverlap));
        }
    }
}

void AfqmcConstraintPath::adjustETAndResetMeasurement()
{
    method.ET = ( observeMeasure.returnEnergy() ).real();

    if( MPIRank()==0 )
    {
        cout<<"\nAdjust trial energy: "<<method.ET<<"\n"<<endl;
    }

    observeMeasure.reSet();
    observeMeasureMixed.reSet();
}

void AfqmcConstraintPath::addMeasurement()
{
    if( !isBP ) return;
    prepareBackPropagationMeasAndNextBackup();

    complex<double> logOverlap;
    double logDiff;
    WalkerRight walkerTemp,walkerTempBP; //WalkerLeft phiLeftOne, phiLeftTwo;
    vector<complex<double>> sumReWeight = returnSumReWeight();

    for(int i = 0; i < method.walkerSizePerThread; ++i)
    {

        if( walkerIsAlive[i] )
        {
            //icf: the ovelap <G|phi'>
            oneBodyWalkerRightOperation.applyToRight(expHalfDtK, walker[i], walkerTemp);   
            walkerWalkerOperation.set(phiT, walkerTemp);   
            logOverlap = walkerWalkerOperation.returnLogOverlap();
            observeMeasureMixed.addMeasurement(walkerWalkerOperation, exp(logOverlap));

            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//icf: phaseless
            //Cap weight real part
            logDiff = ( logOverlap - logOverlapBackup[i] ).real();
            if( (logDiff-method.logEnergyCap)>1e-12 )
            {
                cout<<"Cap weight in measurement: "
                    <<setw(25)<<i+MPIRank()*method.walkerSizePerThread
                    <<setw(25)<<logDiff<<setw(25)<<method.logEnergyCap<<endl;
                logOverlap+= ( method.logEnergyCap - logDiff );
            }

            //Get back all reWeights
            //logOverlap -= sumReWeight[i];   //icf: this sumReWeight is used to add phase at operator BK and should be zero for exact case ?????
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//


            size_t L = walkerTemp.getL();
            TensorHao<std::complex<double>, 2> dotb(L,L);
            TensorHao<std::complex<double>, 2> bforward(L,L);
            TensorHao<std::complex<double>, 2> bbackward(L,L);
            for (size_t j = 1-1; j<=L-1; j++){
                bforward(j,j)=complex <double> (1.0,0.0);
                bbackward(j,j)=complex <double> (1.0,0.0);
            }

            oneBodyWalkerRightOperation.applyToRight(expHalfDtK, walkerBPMeasure[i], walkerTempBP); 
            for (size_t j = 1-1; j<=method.backPropagationStep-1; j++)
            {
                bcsBPMove(i,j,walkerTempBP,bforward,bbackward); 
               
                if( (j)%method.mgsStep == 0 ) walkerTempBP.normalize();
                if( (j)%method.mgsStep == 0 ) stabilize(walkerTempBP,bforward,bbackward,dotb);
            }      

            walkerWalkerOperation.setBCSBP(phiT, walkerTempBP,bforward,bbackward,dotb);

            observeMeasure.addMeasurementBCSBP(walkerWalkerOperation, exp(logOverlap));
        }
    }
    cout<<"Mixed Energy: "<<observeMeasureMixed.returnEnergy()<<endl;
    cout<<"Energy: "<<observeMeasure.returnEnergy()<<endl;
    
}

void AfqmcConstraintPath::bcsBPMove(size_t i, size_t j, WalkerRight &walker, TensorHao<std::complex<double>, 2> &bforward, TensorHao<std::complex<double>, 2> &bbackward)
{
    size_t L = walker.getL();size_t halfL = L/2;
    WalkerRight walkerTemp; 
    TensorHao<std::complex<double>, 2> bforwardTemp(L,L);
    TensorHao<std::complex<double>, 2> bbackwardTemp(L,L);
    double phase1,phase2;
    complex<double> auxForce;

    //T/2
    oneBodyWalkerRightOperation.applyToRight(expMinusHalfDtK, walker, walkerTemp);walker=walkerTemp;
    BL_NAME(gmm)(expMinusHalfDtK.matrix,bforward,bforwardTemp);bforwardTemp *= exp(expMinusHalfDtK.logw);bforward=bforwardTemp;
    BL_NAME(gmm)(expHalfDtK.matrix,bbackward,bbackwardTemp);bbackwardTemp *= exp(expHalfDtK.logw);bbackward=bbackwardTemp;
    walker.addLogw( method.dt*method.ET );

    //V/2
    bforwardTemp.resize(L,L);bbackwardTemp.resize(L,L);
    bforwardTemp=complex <double> (0,0);bbackwardTemp=complex <double> (0,0);

    if( method.forceType == "constForce" ){
        twoBodySample = expMinusDtV.getTwoBodySampleFromAuxForce(twoBodyAuxBPMeasure[i][j], constForce);

        auxForce = 0.0; for(size_t k = 0; k < halfL ; ++k) auxForce += -1.0*twoBodyAuxBPMeasure[i][j](k)*constForce(k);
    }
    else if( method.forceType == "dynamicForce" ){
        walkerWalkerOperation.set(phiT, walker);
        dynamicForce = observeMeasure.getForce(expMinusDtV, walkerWalkerOperation, method.forceCap);
        twoBodySample = expMinusDtV.getTwoBodySampleFromAuxForce(twoBodyAuxBPMeasure[i][j], dynamicForce);

        auxForce = 0.0; for(size_t k = 0; k < halfL ; ++k) auxForce += -1.0*twoBodyAuxBPMeasure[i][j](k)*dynamicForce(k);
        //auxForce=twoBodySample.logw;
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//icf: phaseless
    //Cos projection
    phase2 = auxForce.imag();
    walker.addLogw( log( cos(phase2) ) );
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    twoBodySampleWalkerRightOperation.applyToRight(twoBodySample, walker, walkerTemp);walker=walkerTemp;

    const TensorHao<complex<double>,1> &diag00 = twoBodySample.diag00;
    const TensorHao<complex<double>,1> &diag10 = twoBodySample.diag10;
    const TensorHao<complex<double>,1> &diag01 = twoBodySample.diag01;
    const TensorHao<complex<double>,1> &diag11 = twoBodySample.diag11;

    for(size_t k = 1-1; k <= L-1; k++)
    {
        for(size_t l = 1-1; l <= halfL-1; l++)
        {
            bforwardTemp(l,k)        = diag00(l) * bforward(l,k);// + diag01(l)*bforward(l+halfL, k);
            bforwardTemp(l+halfL, k) = diag11(l) * bforward(l+halfL, k); //diag10(l) * bforward(l,k) + 
        }
    }bforwardTemp *= exp(twoBodySample.logw+log( cos(phase2)) );bforward=bforwardTemp;   //icf: log( cos(phase2)) also apply to bb and bf
    for(size_t k = 1-1; k <= L-1; k++)
    {
        for(size_t l = 1-1; l <= halfL-1; l++)  //icf: ATTENTION: WATCH OUT FOR e^-V !!!!!!!!!!!!!!!!
        {
            bbackwardTemp(l,k)        = conj(1.0/diag00(l)) * bbackward(l,k); // + 1.0/diag10(l)*bbackward(l+halfL, k);
            bbackwardTemp(l+halfL, k) = conj(1.0/diag11(l)) * bbackward(l+halfL, k); //1.0/diag01(l) * bbackward(l,k) +
        }
    }bbackwardTemp *= exp(complex <double> (-1.0,0)*conj(twoBodySample.logw +log(cos(phase2))) );bbackward=bbackwardTemp;
               

    //T/2
    bforwardTemp.resize(L,L);bbackwardTemp.resize(L,L);
    bforwardTemp=complex <double> (0,0);bbackwardTemp=complex <double> (0,0);
    oneBodyWalkerRightOperation.applyToRight(expMinusHalfDtK, walker, walkerTemp);walker=walkerTemp;
    BL_NAME(gmm)(expMinusHalfDtK.matrix,bforward,bforwardTemp);bforwardTemp *= exp(expMinusHalfDtK.logw);bforward=bforwardTemp;
    BL_NAME(gmm)(expHalfDtK.matrix,bbackward,bbackwardTemp);bbackwardTemp *= exp(expHalfDtK.logw);bbackward=bbackwardTemp;  
}

void AfqmcConstraintPath::stabilize(WalkerRight &walker, TensorHao<std::complex<double>, 2> &bforward,TensorHao<std::complex<double>, 2> &bbackward,TensorHao<std::complex<double>, 2> &dotb)
{
    size_t L = walker.getL();size_t N = walker.getN();
    TensorHao<std::complex<double>, 2> projf(N,L);
    TensorHao<std::complex<double>, 2> projb(N,L);

    BL_NAME(gmm)(walker.getWf(),bforward,projf,'C');                      //icf: Attention: we use "C" here instead of "T" (which is used in Ettore's fortran code)
    BL_NAME(gmm)(walker.getWf(),bbackward,projb,'C');

    TensorHao<std::complex<double>, 2> dotbTemp(L,L);
    BL_NAME(gmm)(projb,projf,dotbTemp,'C');
    dotb=dotb+dotbTemp;
    
    TensorHao<std::complex<double>, 2> bforwardNew(L,L);
    TensorHao<std::complex<double>, 2> matrixTemp(L,L);
    TensorHao<std::complex<double>, 2> bbackwardNew(L,L);

    BL_NAME(gmm)(walker.getWf(),projf,matrixTemp);
    bforwardNew=bforward-matrixTemp;

    BL_NAME(gmm)(walker.getWf(),projb,bbackwardNew);

    bforward=bforwardNew;
    bbackward=bbackwardNew;
}


void AfqmcConstraintPath::writeMeasurement()
{
    observeMeasure.write();
    observeMeasureMixed.mixedWrite();
}

void AfqmcConstraintPath::resetMeasurement()
{
    observeMeasure.reSet();
    observeMeasureMixed.reSet();
}
