//
// Created by Hao Shi on 12/18/18.
//

#include "../include/bcs.h"

using namespace std;
using namespace tensor_hao;

void Bcs::readAndBcastModel(const string &filename)
{
    if( MPIRank()==0 )
    {
        ifstream file;
        file.open(filename, ios::in);
        if ( ! file.is_open() ) {cout << "Error opening file in File!!! "<<filename<<endl; exit(1);}

        readFile( L, file );
        readFile( N, file );

        t.resize(L, L); readFile( t.size(),  t.data(),  file );
        D.resize(L, L); readFile( D.size(),  D.data(),  file );
        U.resize(L); readFile( U.size(),  U.data(),  file );

        file.close();

        checkHermitian(t, 1e-6);
        checkHermitian(D, 1e-6);
    }

    MPIBcast( L );
    MPIBcast( N );
    MPIBcast( t );
    MPIBcast( D );
    MPIBcast( U );
}

void Bcs::setH0()
{
    H0.resize(2*L, 2*L);

    for(size_t i = 0; i < L; ++i)
    {
        for(size_t j = 0; j < L; ++j)
        {
            H0(j,   i  ) =  t(j, i)/2.0;
            H0(j+L, i+L) = -conj( t(i, j) )/2.0;

            H0(j+L, i  ) = -D(i,j)/2.0;
            H0(j,   i+L) = -conj( D(j,i) )/2.0;
        }
    }

    checkHermitian(H0, 1e-6);
}

void Bcs::initialHartreeFockEssential()
{
    pairingBefore.resize(L);
    pairing.resize(L);

    if( method.initialType == "setFromModel" )
    {
        pairing = complex<double>(0,0);
    }
    else if( method.initialType == "readWaveFunction" )
    {
        variationalStateF.read("bcsF.dat");
        checkHermitian( variationalStateF );
        calculateVariationalStateFGreen();
        calculatePairing();
        checkPairing();
    }
    else if( method.initialType == "readOrderParameter" )
    {
        pairing.read("pairing.dat");
        checkPairing();
    }
    else
    {
        cout<<"Error!!! Do not know initialType "<<method.initialType<<endl;
        exit(1);
    }

    bcsMu = method.initMu;
    variationalEnergy = 1e300;
}

void Bcs::setMeanFieldVectorsValuesFromPairingAndMu()
{
    //pairing = < Ci_dn^{+} Ci_up^{+} >

    meanFieldVectors = H0;
    for(size_t i = 0; i < L; ++i)
    {
        meanFieldVectors(i,   i  ) += -bcsMu/2.0;
        meanFieldVectors(i+L, i+L) +=  conj(bcsMu)/2.0;

        meanFieldVectors(i+L, i  ) += -U(i) * pairing(i)/2.0;
        meanFieldVectors(i,   i+L) += -U(i) * conj( pairing(i) )/2.0;
    }

    meanFieldValues.resize(2*L);
    BL_NAME(eigen)(meanFieldVectors, meanFieldValues);

//    for(size_t i = 0; i < L; ++i)
//    {
//        for(size_t j = 0; j < L ; ++j)
//        {
//            cout<<"bbb "<<meanFieldVectors(j, i+L)<<endl;
//        }
//    }
}

void Bcs::setVariationalStateFromMeanFieldVectors()
{
    //TODO: CHECK mean field values and X,Y matrix, the other half.
//    cout<<meanFieldValues<<endl;

    TensorHao<std::complex<double>,2> X(L, L), Y(L, L);
    variationalStateF.resize(L, L);

    for(size_t i = 0; i < L; ++i)
    {
        for(size_t j = 0; j < L; ++j)
        {
            X(j, i) = meanFieldVectors(j,   i+L);
            Y(j, i) = meanFieldVectors(j+L, i+L);
        }
    }

    BL_NAME(inverse)(X);
    BL_NAME(gmm)( Y, X, variationalStateF );
    variationalStateF = conj( variationalStateF );

//    cout<<variationalStateF<<endl;
    checkHermitian( variationalStateF, 1e-6 );
//    std::cin.ignore( std::numeric_limits <std::streamsize> ::max(), '\n' );

}

void Bcs::calculateVariationalStateFGreen()
{
    //calculate (1+f^2)^(-1)
    variationalStateFGreen.resize(L, L);
    BL_NAME(gmm)( variationalStateF, variationalStateF, variationalStateFGreen );
    for(size_t i = 0; i < L; ++i)
    {
        variationalStateFGreen(i, i) += 1.0;
    }
    BL_NAME(inverse)(variationalStateFGreen);
}

void Bcs::calculateNparticle()
{
    complex<double> variationalN(0, 0);
    for(size_t i = 0; i < L; ++i)
    {
        variationalN += ( 1.0 - variationalStateFGreen(i, i) );
    }

    if( abs( variationalN.imag() ) > 1e-5 )
    {
        cout<<"Error!!! number of particle of spin up has imaginary part: " << variationalN.imag()<<endl;
        exit(1);
    }

    bcsN = variationalN.real();
}

void Bcs::calculatePairing()
{
    pairingBefore = move( pairing );

    pairing.resize(L);
    for(size_t i = 0; i < L; ++i)
    {
        pairing(i) = 0.0;
        for(size_t j = 0; j < L; ++j)
        {
            pairing(i) += variationalStateFGreen(i, j) * variationalStateF(j, i);
        }
        pairing(i) = conj( pairing(i) );
    }
//    cout<< pairing <<endl;
//    std::cin.ignore( std::numeric_limits <std::streamsize> ::max(), '\n' );
}

void Bcs::calculateVariationalEnergy()
{
    variationalEnergyBefore = variationalEnergy;

    TensorHao<complex<double>,2> Tij(L, L), Dij(L, L);

    Tij = -variationalStateFGreen;
    for(size_t i = 0; i < L; ++i) Tij(i, i) += 1.0;

    BL_NAME(gmm)( variationalStateFGreen, variationalStateF, Dij );
    Dij = conj( Dij );

    variationalEnergy = 0.0;
    for(size_t i = 0; i < L; ++i)
    {
        for(size_t j = 0; j < L; ++j)
        {
            variationalEnergy += 2.0 * ( t(j, i) * Tij(j, i) ).real();  //Both spin, complex conjugate.
            variationalEnergy += 2.0 * ( D(j, i) * Dij(j, i) ).real();  //Pairing part, complex conjugate.
        }
    }
    for(size_t i = 0; i < L; ++i)
    {
        variationalEnergy += U(i) * ( Tij(i,i)*Tij(i,i) + Dij(i,i)*Dij(i,i) ).real();
    }
}
