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
    if( argc != 3 ) {cout<<"Error!!!!! The number of input is not correct"<<endl;}
    string in_filename  = argv[1];
    string out_filename = argv[2];

    vector<complex<double>> vec1;
    size_t Dimen,sizex,sizey;
    double re, im;

    ifstream input_file;
    input_file.open(in_filename, ios::in);
    if ( ! input_file.is_open() ) {cout << "Error opening input file in complex_double_error_analysis.cpp!!!"; exit(1);}
    input_file >> Dimen;
    input_file >> sizex >> sizey;
    while (input_file >> re >> im) vec1.push_back( complex<double>(re, im) );
    input_file.close();

    //write data to output
    ofstream output_file;
    output_file.open(out_filename, ios::out|ios::trunc);
    if( !output_file.is_open() ) {cout << "Error opening output file in complex_double_error_analysis.cpp!!! "<<endl; exit(1);}
    output_file<<setprecision(16)<<scientific;
    for(size_t i=0; i<vec1.size(); i++) output_file<<setw(26)<<real(vec1[i])<<"\n";
    output_file.close();

    return 0;

}
