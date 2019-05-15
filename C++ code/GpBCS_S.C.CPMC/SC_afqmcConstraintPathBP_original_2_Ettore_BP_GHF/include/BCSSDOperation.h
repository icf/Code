//
// Created by boruoshihao on 1/10/17.
// Modefied by icf on 04/29/2019
//

#ifndef AFQMCLAB_BCSSDOPERATION_H
#define AFQMCLAB_BCSSDOPERATION_H

#include "BCS.h"
#include "afqmclab.h"

using namespace std;
using namespace tensor_hao;


enum class BCSSDOperationState
{
    VOID=0,
    BCSSDOVELAPMATRIX,
    THETA_T
};

class BCSSDOperation
{
    BCSSDOperationState state;
    const BCS *walkerLeft;
    const SD *walkerRight;

    tensor_hao::TensorHao< std::complex<double>, 2 > bforward;
    tensor_hao::TensorHao< std::complex<double>, 2 > bbackward;
    tensor_hao::TensorHao< std::complex<double>, 2 > dotB;

    tensor_hao::TensorHao< std::complex<double>, 2 > bcssdOvelapMatrix;       
    tensor_hao::TensorHao< std::complex<double>, 2 > theta_T;      //icf: SD*Q^-1*SD^T

    std::complex<double> logOverlap; bool logOverlapIsCalculated;
    tensor_hao::TensorHao< std::complex<double>, 2 > greenMatrix; bool greenMatrixIsCalculated; bool greenMatrixBCSBPIsCalculated;
    tensor_hao::TensorHao< std::complex<double>, 1 > greenDiagonal; bool greenDiagonalIsCalculated; bool greenDiagonalBCSBPIsCalculated;
    tensor_hao::TensorHao< std::complex<double>, 1 > greenOffDiagonal; bool greenOffDiagonalIsCalculated; bool greenOffDiagonalBCSBPIsCalculated;
    tensor_hao::TensorHao< std::complex<double>, 2 > cmcmMatrix; bool cmcmMatrixIsCalculated; bool cmcmMatrixBCSBPIsCalculated;  //icf: annormal term for BCS <BCS|CC|SD>/<>
    tensor_hao::TensorHao< std::complex<double>, 2 > cpcpMatrix; bool cpcpMatrixIsCalculated; bool cpcpMatrixBCSBPIsCalculated;  //icf: annormal term for BCS <BCS|C^+C^+|SD>/<>
 public:
    BCSSDOperation();
    BCSSDOperation(const BCS &walkerLeft_, const SD &walkerRight_);
    ~BCSSDOperation();

    BCSSDOperationState getState() const;
    const BCS *getWalkerLeft() const;
    const SD *getWalkerRight() const;
    
    void set(const BCS &walkerLeft_, const SD &walkerRight_);
    void setBCSBP(const BCS &walkerLeft_, const SD &walkerRight_, const TensorHao<std::complex<double>, 2> &bforward_, const TensorHao<std::complex<double>, 2> &bbackward_, const TensorHao<std::complex<double>, 2> &dotB_ );
    void reSet();
    
    const tensor_hao::TensorHao< std::complex<double>, 2 > &returnBCSSDOvelapMatrix();
    const tensor_hao::TensorHao<std::complex<double>, 2> &returnTheta_T();
    std::complex<double> returnLogOverlap();
    const tensor_hao::TensorHao< std::complex<double>, 2 > &returnGreenMatrix();
    const tensor_hao::TensorHao< std::complex<double>, 1 > &returnGreenDiagonal();
    const tensor_hao::TensorHao< std::complex<double>, 1 > &returnGreenOffDiagonal();
    const tensor_hao::TensorHao< std::complex<double>, 2 > &returnCmCmMatrix();
    const tensor_hao::TensorHao< std::complex<double>, 2 > &returnCpCpMatrix();

    const tensor_hao::TensorHao< std::complex<double>, 2 > &returnGreenMatrixBCSBP();
    const tensor_hao::TensorHao< std::complex<double>, 1 > &returnGreenDiagonalBCSBP();
    const tensor_hao::TensorHao< std::complex<double>, 1 > &returnGreenOffDiagonalBCSBP();
    const tensor_hao::TensorHao< std::complex<double>, 2 > &returnCmCmMatrixBCSBP();
    const tensor_hao::TensorHao< std::complex<double>, 2 > &returnCpCpMatrixBCSBP();

    double getMemory() const;

 private:
    BCSSDOperation(const BCSSDOperation& x);
    BCSSDOperation & operator  = (const BCSSDOperation& x);

    void calculateBCSSDOvelapMatrix();
    void calculateTheta_T();
};

void setWalkerFromPhiT(SD &walker, const BCS& phiT, size_t L_temp, size_t N_temp);     

#endif //AFQMCLAB_BCSSDOPERATION_H
