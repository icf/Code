#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <complex>
#include <typeinfo>
#include <tuple>

using namespace std;

int main(int argc, char** argv)
{
    if( argc != 4 ) {cout<<"Error!!!!! The number of input is not coorect"<<endl;}
    string in_filename  = argv[1];
    string out_filename = argv[2];    //spin
    string out_filename2 = argv[3];   //charge

    vector<complex<double>> vec;
    vector<complex<double>> spinDensity;
    vector<complex<double>> chargeDensity;

    //read data to vec
    double re, im;
    int dimen, Lx, Ly;

    ifstream input_file;
    input_file.open(in_filename, ios::in);
    if ( ! input_file.is_open() ) {cout << "Error opening input file in complex_double_error_analysis.cpp!!!"; exit(1);}
    
    input_file >> dimen;input_file >> Lx;input_file >> Ly;
    while (input_file >> re >> im) vec.push_back( complex<double>(re, im) );
    input_file.close();


    cout<<"Green Matrix is "<<vec.size()<<" with standard: "<<Lx*Ly<<endl;

    //data analysis
    
    

    //write data to output
    ofstream output_file;
    output_file.open(out_filename, ios::out|ios::trunc);
    if( !output_file.is_open() ) {cout << "Error opening output file in complex_double_error_analysis.cpp!!! "<<endl; exit(1);}
    output_file<<setprecision(16)<<scientific;
    for(size_t i=1-1; i<=Lx/2-1; i++) output_file<<setw(26)<< real(vec[i+i*Lx]-vec[i+Lx/2+(i+Lx/2)*Lx])/2.0 <<"\n";
    output_file.close();

    ofstream output_file2;
    output_file2.open(out_filename2, ios::out|ios::trunc);
    if( !output_file2.is_open() ) {cout << "Error opening output file in complex_double_error_analysis.cpp!!! "<<endl; exit(1);}
    output_file2<<setprecision(16)<<scientific;
    for(size_t i=1-1; i<=Lx/2-1; i++) output_file2<<setw(26)<< real(vec[i+i*Lx]+vec[i+Lx/2+(i+Lx/2)*Lx])/2.0 <<"\n";
    output_file2.close();

    return 0;

}
