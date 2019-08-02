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

    string outputname = "densityMatrixTensorHao.dat";

    size_t Dimen,sizex,sizey;

    //read green matrix
    ifstream input_file;
    input_file.open(in_filename, ios::in);
    if ( ! input_file.is_open() ) {cout << "Error opening input file in complex_double_error_analysis.cpp!!!"<<endl; exit(1);}
    input_file >> Dimen;
    input_file >> sizex >> sizey;

    //MatrixXcd cicj_global=MatrixXcd::Zero(sizex,sizey);
    TensorHao< complex<double>, 2 > cicj_global(sizex,sizey);cicj_global=complex<double>(0.0,0.0);
    for(size_t i=1-1;i<=sizey-1;i++){
       for(size_t j=1-1;j<=sizex-1;j++){
          input_file >> cicj_global(j,i);       
       }
    }  
    input_file.close();

    TensorHao< complex<double>, 2 > densityMatrix(sizex,sizey);
    densityMatrix=complex <double> (0.0,0.0);

    for(size_t i=1-1;i<=sizey-i;i++){
       for(size_t j=1-1;j<=sizex-1;j++){
          densityMatrix(j,i)=cicj_global(j,i);
       }
    } 

    //write data to output
    ofstream file;
    file.open(outputname, ios::out|ios::trunc);
    if ( ! file.is_open() ) { cout << "Error opening file in File!!! "<<outputname<<endl; exit(1); }
    densityMatrix.write(file);
    file.close();

    return 0;

}
