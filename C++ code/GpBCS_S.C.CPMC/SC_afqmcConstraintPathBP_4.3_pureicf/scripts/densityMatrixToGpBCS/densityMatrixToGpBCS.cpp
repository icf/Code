#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <complex>
#include <typeinfo>
#include <tuple>
#include "afqmclab.h"

using namespace std;
using namespace tensor_hao;

using namespace std;

int main(int argc, char** argv)
{
    if( argc != 2 ) {cout<<"Error!!!!! The number of input is not correct"<<endl;}
    string in_filename  = argv[1];

    string wf = "pBCS_PhiT.Wf.dat";
    string T = "pBCS_PhiT.T.dat";
    string Gdia = "pBCS_PhiT.Gdia.dat";

    size_t Dimen,sizex,sizey,L;
    double re, im;

    //read green matrix
    ifstream input_file;
    input_file.open(in_filename, ios::in);
    if ( ! input_file.is_open() ) {cout << "Error opening input file in complex_double_error_analysis.cpp!!!"<<endl; exit(1);}
    input_file >> Dimen;
    input_file >> sizex >> sizey;

    if(sizex == sizey){
      L=sizex;
    }else{
      cout << "Error green matrix size!!!"<<endl;; 
      exit(1);
    }

    TensorHao< complex<double>, 2 > greenMatrix(L,L),tempOrbital(L,L),tempAnalytic(L,L);
    greenMatrix=complex <double> (0.0,0.0);tempOrbital=complex <double> (0.0,0.0);tempAnalytic=complex <double> (0.0,0.0);
    TensorHao< double, 1 > tempDia(L);tempDia=0.0;


    for(size_t j=1-1;j<=L-1;j++){
       for(size_t i=1-1;i<=L-1;i++){
          input_file >> re >> im;
          greenMatrix(i,j)=( complex<double>(re, im) );
       }
    } 
    input_file.close();

    //make phiT
    checkHermitian(greenMatrix); 

    TensorHao<complex<double>, 2> testMatrix(L,L),testMatrix2(L,L),greenMatrixSave(L,L);
    greenMatrixSave=greenMatrix;

    BL_NAME(eigen)(greenMatrix,tempDia);
    tempOrbital=greenMatrix;
    cout<<greenMatrix<<endl;
 
    cout<<tempDia<<endl;
    BL_NAME(gmm)(tempOrbital, greenMatrixSave, testMatrix, 'C'); 
    BL_NAME(gmm)(testMatrix, tempOrbital, testMatrix2); 
    cout<<testMatrix2<<endl;exit(1);

    cout<<"eigenValue: "<<tempDia<<endl;
   
    for(size_t i=1-1; i<=L/2-1; i++){
        complex <double> temp=complex <double>((tempDia(2*i)+tempDia(2*i+1))/2.0);
        tempAnalytic(2*i,2*i+1)=sqrt(temp/(complex <double>(1,0)-temp));
        tempAnalytic(2*i+1,2*i)=complex<double>(-1.0,0.0)*tempAnalytic(2*i,2*i+1);
    }

    TensorHao<complex<double>, 2> tempMatrix(L,L);
    TensorHao<complex<double>, 2> tempMatrix2(L,L);

    BL_NAME(gmm)(tempOrbital, tempAnalytic, tempMatrix ); 
    BL_NAME(gmm)(tempMatrix, tempOrbital, tempMatrix2 ,'N','T'); 

    TensorHao< complex<double>, 2 > phiT_wfRef=tempMatrix2;
    complex<double> phiT_logwRef=complex<double>(0.0,0.0);
    TensorHao< complex<double>, 2 > phiT_orbital=tempOrbital;
    TensorHao< complex<double>, 2 > phiT_occupancy=tempAnalytic;

    //write data to output
    ofstream file1;
    file1.open(wf, ios::out|ios::trunc);
    if ( ! file1.is_open() ) { cout << "Error opening file in File!!! "<<wf<<endl; exit(1); }
    writeFile(phiT_logwRef, file1);
    phiT_wfRef.write(file1);
    file1.close();

    ofstream file2;
    file2.open(T, ios::out|ios::trunc);
    if ( ! file2.is_open() ) { cout << "Error opening file in File!!! "<<T<<endl; exit(1); }
    writeFile(phiT_logwRef, file2);
    phiT_orbital.write(file2);
    file2.close();

    ofstream file3;
    file3.open(Gdia, ios::out|ios::trunc);
    if ( ! file3.is_open() ) { cout << "Error opening file in File!!! "<<Gdia<<endl; exit(1); }
    writeFile(phiT_logwRef, file3);
    phiT_occupancy.write(file3);
    file3.close();



    return 0;

}
