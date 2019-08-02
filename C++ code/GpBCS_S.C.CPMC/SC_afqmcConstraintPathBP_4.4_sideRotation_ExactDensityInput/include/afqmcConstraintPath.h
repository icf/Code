//
// Created by boruoshihao on 4/16/17.
// Modefied by icf on 05/02/2019
//

#ifndef AFQMCLAB_AFQMCCONSTRAINTPATH_H
#define AFQMCLAB_AFQMCCONSTRAINTPATH_H

#include "afqmcConstraintPathDefine.h"
#include "afqmcConstraintPathMethod.h"

using namespace std;
using namespace tensor_hao;

class AfqmcConstraintPath
{
 private:
    AfqmcConstraintPathMethod method;
    Model model;
    OneBody expMinusDtK, expMinusHalfDtK, expHalfDtK;
    TwoBody expMinusDtV;
    TwoBodyForce dynamicForce, constForce;

    OneBodyWalkerRightOperation oneBodyWalkerRightOperation;
    OneBodyWalkerLeftOperation oneBodyWalkerLeftOperation;
    TwoBodySampleWalkerRightOperation twoBodySampleWalkerRightOperation;
    TwoBodySampleWalkerLeftOperation twoBodySampleWalkerLeftOperation;

    TwoBodyAux twoBodyAux;
    TwoBodySample twoBodySample;

    WalkerLeft phiT;
    std::vector<WalkerRight> walker;
    std::vector<bool> walkerIsAlive;

    WalkerWalkerOperation walkerWalkerOperation;
    ModelMeasureObserve observeMeasure;
    ModelMeasureObserve observeMeasureMixed;   //only used in BP measure (not mixed ET in thermal)

    bool isBP;
    size_t BPIndex;
    std::vector<WalkerRight> walkerBackup, walkerBPMeasure;
    std::vector<std::vector<TwoBodyAux>> twoBodyAuxBackup, twoBodyAuxBPMeasure;
    std::vector<std::vector<int>> tablePerThreadBackup;
 public:
    AfqmcConstraintPath();
    ~AfqmcConstraintPath();

    void run();
    void runSC();
    void initialParameters();
    void initialMeasure();
    void estimateMemory();
    void measureWithoutProjection();
    void measureWithProjection();
    void prepareStop();
    // for BCS BK, original from Ettore's fortran code
    void bcsBPMove(size_t i, size_t j, WalkerRight &walker, TensorHao<std::complex<double>, 2> &bforward, TensorHao<std::complex<double>, 2> &bbackward);
    void stabilize(WalkerRight &walker, TensorHao<std::complex<double>, 2> &bforward,TensorHao<std::complex<double>, 2> &bbackward,TensorHao<std::complex<double>, 2> &dotB);
    
    //cut line flag
    size_t scLooPSteps; // the S.C. step.
    bool cutMinus;

    //icf: get the energy and green matrix after one CPMC 
    std::vector<TensorHao<std::complex<double>, 2>> greenMatrixAfterCPMCMeasureSave;
    TensorHao<std::complex<double>, 2> greenMatrixAfterCPMC,greenMatrixAfterCPMCErrorBar;  
    std::vector<TensorHao<std::complex<double>, 2>> densityDensityMatrixAfterCPMCMeasureSave;
    TensorHao<std::complex<double>, 2> densityDensityMatrixAfterCPMC,densityDensityMatrixAfterCPMCErrorBar;  
    void saveMeasureGreenMatrix();
    void resetAndSaveGreenMatrix();

 private:
    void initialPhiT();
    void initialSCPhiT();
    void initialWalker();
    void writeWalkers();
    void initialMgsAndPopControl();

    void projectExpHalfDtK();
    void projectExpMinusHalfDtK();
    void projectExpMinusDtKExpMinusDtV();
    void projectOneStep(size_t &mgsIndex, size_t &popControlIndex);
    void modifyGM();
    void popControl();
    void checkAndResetWalkerIsAlive();

    void addMixedMeasurement();
    void adjustETAndResetMeasurement();
    void addMeasurement();
    void writeMeasurement();
    void resetMeasurement();

    void prepareBackPropagation();
    void prepareBackPropagationMeasAndNextBackup();
    void generateWalkerBackupFromWalkerBPMeasure();
};

#endif //AFQMCLAB_AFQMCCONSTRAINTPATH_H
