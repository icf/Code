//
// Created by Hao Shi on 12/18/18.
//

#ifndef BCSHUBBARD_BCS_H
#define BCSHUBBARD_BCS_H

#include "afqmclab.h"
#include "bcsMethod.h"

class Bcs
{
 private:
    BcsMethod method;

    size_t L, N;
    tensor_hao::TensorHao<std::complex<double>, 2> t, D;
    tensor_hao::TensorHao<double, 1> U;
    tensor_hao::TensorHao<std::complex<double>, 2> H0;

    double variationalEnergy, variationalEnergyBefore;
    tensor_hao::TensorHao<std::complex<double>,2> variationalStateF, variationalStateFGreen;
    tensor_hao::TensorHao<std::complex<double>,1> pairing, pairingBefore;
    tensor_hao::TensorHao<std::complex<double>, 2> meanFieldVectors;
    tensor_hao::TensorHao<double, 1> meanFieldValues;

    double bcsMu, bcsMuPlus, bcsMuMinus;
    double bcsN;

    double convergeValue;
    double minimumEnergy;
    tensor_hao::TensorHao<std::complex<double>,2>  minimumStateF;

 public:
    Bcs();
    ~Bcs();

    void run();
    void initialParameters();
    void prepareStop();

 private:
    void readAndBcastModel(const std::string &filename);
    void setH0();
    void initialHartreeFockEssential();
    void setMeanFieldVectorsValuesFromPairingAndMu();
    void setVariationalStateFromMeanFieldVectors();
    void calculateVariationalStateFGreen();
    void calculateNparticle();
    void calculatePairing();
    void calculateVariationalEnergy();

    void findPairing();
    void calculateConvergeValue();
    void relaxPairing();
    void annealPairing(double anneal);
    void checkPairing();

    void findBcsMu();
    void calculateBcsNFromPairingAndMu();
    void findBcsMuPlusBcsMuMinus();
    void findBcsMuBetweenPlusAndMinus();
};
#endif //BCSHUBBARD_BCS_H
