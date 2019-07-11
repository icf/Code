//
// Created by boruoshihao on 1/10/17.
// Modefied by icf on 04/29/2019
//

#include "../include/BCSSDOperation.h"
#include "../include/afqmcConstraintPath.h"


using namespace std;
using namespace tensor_hao;

BCSSDOperation::BCSSDOperation() : walkerLeft(nullptr), walkerRight(nullptr)
{
    reSet();
}

BCSSDOperation::BCSSDOperation(const BCS &walkerLeft_, const SD &walkerRight_)
{
    set(walkerLeft_, walkerRight_);
}

BCSSDOperation::~BCSSDOperation() { }

BCSSDOperationState BCSSDOperation::getState() const { return state; }

const BCS *BCSSDOperation::getWalkerLeft() const { return walkerLeft; }

const SD *BCSSDOperation::getWalkerRight() const { return walkerRight; }

void BCSSDOperation::set(const BCS &walkerLeft_, const SD &walkerRight_)
{
    walkerLeft  = &walkerLeft_;
    walkerRight = &walkerRight_;
    size_t L = walkerLeft->getL(); 
    if( L != walkerRight->getL() )
    {
        cout<<"Error!!! Find Inconsistency between walkerLeft and walkerRight!"<<endl;
        exit(1);
    }

    reSet();
}

void BCSSDOperation::setBCSBP(const BCS &walkerLeft_, const SD &walkerRight_, const TensorHao<std::complex<double>, 2> &bforward_, const TensorHao<std::complex<double>, 2> &bbackward_, const TensorHao<std::complex<double>, 2> &dotB_ )
{
    walkerLeft  = &walkerLeft_;
    walkerRight = &walkerRight_;

    bforward=bforward_;
    bbackward=bbackward_;
    dotB=dotB_;

    size_t L = walkerLeft->getL(); 
    if( L != walkerRight->getL() )
    {
        cout<<"Error!!! Find Inconsistency between walkerLeft and walkerRight!"<<endl;
        exit(1);
    }

    reSet();
}

void BCSSDOperation::reSet()
{
    state = BCSSDOperationState::VOID;
    logOverlapIsCalculated = false;
    greenMatrixIsCalculated = false;
    greenDiagonalIsCalculated = false;
    greenOffDiagonalIsCalculated = false;
    cmcmMatrixIsCalculated = false;
    cpcpMatrixIsCalculated = false;
    greenMatrixBCSBPIsCalculated = false;
    greenDiagonalBCSBPIsCalculated = false;
    greenOffDiagonalBCSBPIsCalculated = false;
    cmcmMatrixBCSBPIsCalculated = false;
    cpcpMatrixBCSBPIsCalculated = false;
}

const TensorHao<complex<double>, 2> &BCSSDOperation::returnBCSSDOvelapMatrix()
{
    calculateBCSSDOvelapMatrix();
    return bcssdOvelapMatrix;
}

const TensorHao<complex<double>, 2> &BCSSDOperation::returnTheta_T()   ////icf: Theta_T-->SD*Q^-1*SD^T
{
    calculateTheta_T();
    return theta_T;
}

complex<double> BCSSDOperation::returnLogOverlap()
{
    if( logOverlapIsCalculated ) return logOverlap;
    calculateBCSSDOvelapMatrix();
    checkSkewSymmetric(bcssdOvelapMatrix,1e-10);
    TensorHao<complex<double>, 2> tempMatrix=bcssdOvelapMatrix;
    logOverlap = conj(walkerLeft->getLogw()) + walkerRight->getLogw() + log(pfaffian(tempMatrix));
    logOverlapIsCalculated = true;
    return logOverlap;
}

const TensorHao<complex<double>, 2> &BCSSDOperation::returnGreenMatrix()
{
    if(greenMatrixIsCalculated) return greenMatrix;
    
    calculateTheta_T();

    size_t L = walkerLeft->getL();
    greenMatrix.resize(L,L);greenMatrix=complex <double> (0,0);

    BL_NAME(gmm)(walkerLeft->getWf(), theta_T , greenMatrix ,'C','N'); 

    greenMatrixIsCalculated = true;
    
    return greenMatrix;
}

const TensorHao<complex<double>, 1> &BCSSDOperation::returnGreenDiagonal()
{
    if(greenDiagonalIsCalculated) return greenDiagonal;
    
    calculateTheta_T();

    size_t L = walkerLeft->getL(); 
    const TensorHao<complex<double>, 2> &wfLeft = walkerLeft->getWf();

    greenDiagonal.resize(L); greenDiagonal = complex<double>(0,0);
    for(size_t j = 0; j < L; ++j)
    {
        for(size_t i = 0; i < L; ++i)
        {
            greenDiagonal(i) += theta_T(j,i) * conj( wfLeft(j, i) );
        }
    }

    greenDiagonalIsCalculated = true;
    
    return greenDiagonal;
}

const TensorHao<complex<double>, 1> &BCSSDOperation::returnGreenOffDiagonal()
{
    if(greenOffDiagonalIsCalculated) return greenOffDiagonal;    

    calculateTheta_T();

    size_t L = walkerLeft->getL();
    const TensorHao<complex<double>, 2> &wfLeft = walkerLeft->getWf();

    size_t halfL = L/2;
    if( L != halfL*2 ) { cout<<"Error!!! Green Matrix rank size is odd number! "<<L<<endl; exit(1); }

    greenOffDiagonal.resize(L); greenOffDiagonal = complex<double>(0,0);

    for(size_t j = 0; j < L; ++j)
    {
        for(size_t i = 0; i < halfL; ++i)
        {
            greenOffDiagonal(i)       += theta_T(j,i) * conj( wfLeft(j, i+halfL) ); // C_i^\dagger C_{i+halfL}
            greenOffDiagonal(i+halfL) += theta_T(j,i+halfL) * conj( wfLeft(j, i) ); // C_{i+halfL}^\dagger C_i
        }
    }
    
    greenOffDiagonalIsCalculated = true;
    
    return greenOffDiagonal;
}

const TensorHao<complex<double>, 2> &BCSSDOperation::returnCmCmMatrix()
{
    if(cmcmMatrixIsCalculated) return cmcmMatrix;
    
    calculateTheta_T();

    size_t L = walkerLeft->getL();
    cmcmMatrix.resize(L,L);cmcmMatrix=complex <double> (0,0);

    cmcmMatrix=complex <double> (-1.0,0)*theta_T;    

    cmcmMatrixIsCalculated = true;
    
    return cmcmMatrix;
}

const TensorHao<complex<double>, 2> &BCSSDOperation::returnCpCpMatrix()
{
    if(cpcpMatrixIsCalculated) return cpcpMatrix;
    
    calculateTheta_T();

    size_t L = walkerLeft->getL();
    cpcpMatrix.resize(L,L);cpcpMatrix=complex <double> (0,0);
    TensorHao<complex<double>,2> tempMatrix(L,L);

    BL_NAME(gmm)(walkerLeft->getWf(), theta_T, tempMatrix ,'C');  
    BL_NAME(gmm)(tempMatrix, walkerLeft->getWf(), cpcpMatrix ,'N','C');  

    cpcpMatrix=complex <double> (-1.0,0)*conjtrans(walkerLeft->getWf()) + cpcpMatrix;    

    cpcpMatrixIsCalculated = true;
    
    return cpcpMatrix;
}

const TensorHao<complex<double>, 2> &BCSSDOperation::returnGreenMatrixBCSBP()     //icf: the formula for Green function may has a Trans-error (we haven't correct it for now at GreenDia and GreenOffDia)
{
    if(greenMatrixBCSBPIsCalculated) return greenMatrix;
    
    calculateTheta_T();

    size_t L = walkerLeft->getL();
    greenMatrix.resize(L,L);greenMatrix=complex <double> (0,0);
    TensorHao< complex<double>, 2 > greenMatrixOne(L,L);greenMatrixOne=complex <double> (0,0);
    TensorHao< complex<double>, 2 > greenMatrixTwo(L,L);greenMatrixTwo=complex <double> (0,0);

    BL_NAME(gmm)(walkerLeft->getWf(), theta_T  , greenMatrixOne ,'C'); 
    BL_NAME(gmm)(bforward, greenMatrixOne , greenMatrixTwo, 'T' ); 
    BL_NAME(gmm)(greenMatrixTwo, conj( bbackward ) , greenMatrix); 

    greenMatrix += trans(dotB);               //icf: ATTENTION for stabilization

    greenMatrixBCSBPIsCalculated = true;
    
    return greenMatrix;
}

const TensorHao<complex<double>, 1> &BCSSDOperation::returnGreenDiagonalBCSBP()
{
    if(greenDiagonalBCSBPIsCalculated) return greenDiagonal;

    calculateTheta_T();

    size_t L = walkerLeft->getL();
    greenMatrix.resize(L,L);greenMatrix=complex <double> (0,0);
    TensorHao< complex<double>, 2 > greenMatrixOne(L,L);
    TensorHao< complex<double>, 2 > greenMatrixTwo(L,L);

    BL_NAME(gmm)(theta_T, conj( walkerLeft->getWf() ) , greenMatrixOne ,'T');  
    BL_NAME(gmm)(bforward, greenMatrixOne , greenMatrixTwo,'T' ); 
    BL_NAME(gmm)(greenMatrixTwo, conj( bbackward ) , greenMatrix); 

    greenMatrix += trans(dotB);

    greenDiagonal.resize(L);greenDiagonal = complex<double>(0,0);
    for(size_t i=1-1;i<=L-1;i++)greenDiagonal(i)=greenMatrix(i,i);

    greenMatrixBCSBPIsCalculated = true;
    
    greenDiagonalBCSBPIsCalculated = true;
    
    return greenDiagonal;
}

const TensorHao<complex<double>, 1> &BCSSDOperation::returnGreenOffDiagonalBCSBP()
{
    if(greenOffDiagonalBCSBPIsCalculated) return greenOffDiagonal;    

    calculateTheta_T();

    size_t L = walkerLeft->getL();
    size_t halfL=L/2;
    greenMatrix.resize(L,L);greenMatrix=complex <double> (0,0);
    TensorHao< complex<double>, 2 > greenMatrixOne(L,L);
    TensorHao< complex<double>, 2 > greenMatrixTwo(L,L);

    BL_NAME(gmm)(theta_T, conj( walkerLeft->getWf() ) , greenMatrixOne ,'T');  
    BL_NAME(gmm)(bforward, greenMatrixOne , greenMatrixTwo,'T' ); 
    BL_NAME(gmm)(greenMatrixTwo, conj( bbackward ) , greenMatrix); 

    greenMatrix += trans(dotB);
    
    greenOffDiagonal.resize(L); greenOffDiagonal = complex<double>(0,0);
    for(size_t i = 0; i < halfL; ++i)
    {
        greenOffDiagonal(i) = greenMatrix(i,i+halfL);
        greenOffDiagonal(i+halfL) = greenMatrix(i+halfL,i);
    }

    greenOffDiagonalBCSBPIsCalculated = true;
    
    return greenOffDiagonal;
}

const TensorHao<complex<double>, 2> &BCSSDOperation::returnCmCmMatrixBCSBP()
{
    if(cmcmMatrixBCSBPIsCalculated) return cmcmMatrix;
    
    calculateTheta_T();

    size_t L = walkerLeft->getL();
    cmcmMatrix.resize(L,L);cmcmMatrix=complex <double> (0,0);
    TensorHao< complex<double>, 2 > cmcmMatrixOne(L,L);
    TensorHao< complex<double>, 2 > cmcmMatrixTwo(L,L);

    cmcmMatrixOne=complex <double> (-1.0,0)*theta_T;      
    BL_NAME(gmm)(bbackward, cmcmMatrixOne , cmcmMatrixTwo,'C' ); 
    BL_NAME(gmm)(cmcmMatrixTwo, conj( bbackward ) , cmcmMatrix); 

    cmcmMatrixBCSBPIsCalculated = true;
    
    return cmcmMatrix;
}

const TensorHao<complex<double>, 2> &BCSSDOperation::returnCpCpMatrixBCSBP()
{
    if(cpcpMatrixBCSBPIsCalculated) return cpcpMatrix;
    
    calculateTheta_T();

    size_t L = walkerLeft->getL();
    cpcpMatrix.resize(L,L);cpcpMatrix=complex <double> (0,0);
    TensorHao<complex<double>,2> tempMatrix(L,L);
    TensorHao< complex<double>, 2 > cpcpMatrixOne(L,L);
    TensorHao< complex<double>, 2 > cpcpMatrixTwo(L,L);

    BL_NAME(gmm)(walkerLeft->getWf(), theta_T, tempMatrix ,'C');  
    BL_NAME(gmm)(tempMatrix, walkerLeft->getWf(), cpcpMatrixOne ,'N','C');   
    cpcpMatrixOne=complex <double> (-1.0,0)*conjtrans(walkerLeft->getWf()) + cpcpMatrixOne;   

    BL_NAME(gmm)(bforward, cpcpMatrixOne , cpcpMatrixTwo,'T' ); 
    BL_NAME(gmm)(cpcpMatrixTwo, bforward  , cpcpMatrix);  

    cpcpMatrixBCSBPIsCalculated = true;
    
    return cpcpMatrix;
}

double BCSSDOperation::getMemory() const
{
    return 8.0*2+bcssdOvelapMatrix.getMemory()+theta_T.getMemory()
           +16.0+1.0+greenMatrix.getMemory()+1.0 +greenDiagonal.getMemory()+1.0
           +greenOffDiagonal.getMemory()+1.0;
}

BCSSDOperation::BCSSDOperation(const BCSSDOperation &x) { }

BCSSDOperation &BCSSDOperation::operator=(const BCSSDOperation &x) { return *this; }

void BCSSDOperation::calculateBCSSDOvelapMatrix()
{
    if( state >= BCSSDOperationState::BCSSDOVELAPMATRIX ) return;

    size_t L = walkerLeft->getL();size_t N = walkerRight->getN();
    TensorHao<complex<double>,2> overlapMatrix(N,N);
    TensorHao<complex<double>,2> tempMatrix(N,L);
    BL_NAME(gmm)( walkerRight->getWf(), conj(walkerLeft->getWf()),  tempMatrix ,'T');
    BL_NAME(gmm)( tempMatrix, walkerRight->getWf(), overlapMatrix );
    bcssdOvelapMatrix = move(overlapMatrix);

    state = BCSSDOperationState::BCSSDOVELAPMATRIX;
}

void BCSSDOperation::calculateTheta_T()
{
    if( state >= BCSSDOperationState::THETA_T ) return;

    size_t L = walkerLeft->getL();size_t N = walkerRight->getN();
    TensorHao<complex<double>,2> qMatrix(N,N);
    TensorHao<complex<double>,2> tempMatrix(N,L);
    TensorHao<complex<double>,2> tempMatrixT(L,N);
    theta_T.resize(L,L);theta_T=complex <double> (0,0);

    //qMatrix=SD^T*BCS^\dagger*SD
    BL_NAME(gmm)( walkerRight->getWf(), walkerLeft->getWf(),  tempMatrix ,'T','C');
    BL_NAME(gmm)( tempMatrix, walkerRight->getWf(),  qMatrix );

    //Theta_T=SD*qMatrix^-1*SD^T
    LUDecomp< std::complex<double> > LU_qMatrix=BL_NAME(LUconstruct)(qMatrix);
    TensorHao<complex<double>,2> inverse_qMatrix=BL_NAME(inverse)(LU_qMatrix);

    BL_NAME(gmm)( walkerRight->getWf(), inverse_qMatrix,  tempMatrixT );
    BL_NAME(gmm)( tempMatrixT, walkerRight->getWf(),  theta_T,'N','T' );

    state = BCSSDOperationState::THETA_T;
}

void setWalkerFromPhiT(SD &walker, const BCS& phiT, size_t L, size_t N)
{

     TensorHao<complex<double>,2> tempMatrix(L,N);tempMatrix=complex <double> (0,0);
     for(size_t k=L-N;k<=L-1;k++){
        cout<<"we take: "<<phiT.occupancy(k,k-1)<<" or "<<phiT.occupancy(k-1,k)<<" as 1"<<endl;
     for(size_t j=1-1;j<=L-1;j++){
        tempMatrix(j,k-L+N) = phiT.orbital(j,k);
     }
     }

     walker.wfRef()=tempMatrix;
     walker.logwRef()=complex <double> (0,0);
}
