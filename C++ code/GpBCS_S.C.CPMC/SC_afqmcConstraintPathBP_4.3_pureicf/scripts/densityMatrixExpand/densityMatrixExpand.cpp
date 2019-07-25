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

    string outputname = "expandedDensityMatrix.dat";

    size_t Dimen,sizex,sizey,L;
    double re, im;

    //read green matrix
    ifstream input_file;
    input_file.open(in_filename, ios::in);
    if ( ! input_file.is_open() ) {cout << "Error opening input file in complex_double_error_analysis.cpp!!!"<<endl; exit(1);}
    input_file >> Dimen;
    input_file >> sizex >> sizey;

    if(sizex == sizey){
      L=2*sizex;
    }else{
      cout << "Error green matrix size!!!"<<endl;; 
      exit(1);
    }

    TensorHao< complex<double>, 2 > densityMatrix(L/2,L/2),expandedDensityMatrix(L,L);
    densityMatrix=complex <double> (0.0,0.0);expandedDensityMatrix=complex <double> (0.0,0.0);

    for(size_t j=1-1;j<=L/2-1;j++){
       for(size_t i=1-1;i<=L/2-1;i++){
          input_file >> re >> im;
          densityMatrix(i,j)=( complex<double>(re, im) );
       }
    } 
    input_file.close();

    //make phiT
    checkHermitian(densityMatrix);

    for(size_t j=1-1;j<=L/2-1;j++){
       for(size_t i=1-1;i<=L/2-1;i++){
          expandedDensityMatrix(i,j)=densityMatrix(i,j);
          expandedDensityMatrix(i+L/2,j+L/2)=expandedDensityMatrix(i,j);
       }
    } 


    //write data to output
    ofstream file;
    file.open(outputname, ios::out|ios::trunc);
    if ( ! file.is_open() ) { cout << "Error opening file in File!!! "<<outputname<<endl; exit(1); }
    expandedDensityMatrix.write(file);
    file.close();




    return 0;

}
