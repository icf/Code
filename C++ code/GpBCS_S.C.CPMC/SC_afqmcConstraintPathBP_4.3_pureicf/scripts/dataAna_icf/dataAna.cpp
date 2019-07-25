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
    string in_filename1  = argv[1];
    string in_filename2  = argv[2];
    string out_filename = argv[3];
    size_t anaStep = atoi(argv[4]);

    vector<complex<double>> vec1,vec2;
    vector<size_t> fact;
    vector<complex<double>> mean,errorBar;

    //read data to vec
    double re, im;

    ifstream input_file;
    input_file.open(in_filename1, ios::in);
    if ( ! input_file.is_open() ) {cout << "Error opening input file in complex_double_error_analysis.cpp!!!"; exit(1);}
    while (input_file >> re >> im) vec1.push_back( complex<double>(re, im) );
    input_file.close();

    input_file.open(in_filename2, ios::in);
    if ( ! input_file.is_open() ) {cout << "Error opening input file in complex_double_error_analysis.cpp!!!"; exit(1);}
    while (input_file >> re >> im) vec2.push_back( complex<double>(re, im) );
    input_file.close();

    cout<<"Effective sample points is "<<vec1.size()<<" and "<<vec2.size()<<endl;
    if(vec1.size() != vec2.size())cout<<"Error!!!! vec1.size() != vec2.size()"<<endl;

    //data analysis
    for(size_t i=1-1; i<=size_t(vec1.size()/anaStep) -1 ;i++ ){
       fact.push_back( (i+1)*anaStep );
       complex<double> sumH=complex<double>(0.0,0.0);complex<double> sumD=complex<double>(0.0,0.0);
       complex<double> tempH=complex<double>(0.0,0.0);complex<double> tempD=complex<double>(0.0,0.0);complex<double> sumErr=complex<double>(0.0,0.0);
       for(size_t j=1-1; j<=anaStep -1; j++ ){
           sumH += vec1[i*anaStep + j];
           sumD += vec2[i*anaStep + j];    
       }
       for(size_t j=1-1; j<=anaStep -1; j++ ){
           tempH = vec1[i*anaStep + j];
           tempD = vec2[i*anaStep + j];    
           sumErr += (tempH/tempD - sumH/sumD)*(tempH/tempD - sumH/sumD);
       }
       mean.push_back(sumH/sumD);
       errorBar.push_back(sqrt(sumErr)/complex <double>(anaStep));
    }

    //write data to output
    ofstream output_file;
    output_file.open(out_filename, ios::out|ios::trunc);
    if( !output_file.is_open() ) {cout << "Error opening output file in complex_double_error_analysis.cpp!!! "<<endl; exit(1);}
    output_file<<setprecision(16)<<scientific;
    for(size_t i=0; i<fact.size(); i++) output_file<<setw(26)<<fact[i] <<setw(26)<<mean[i].real()<<setw(26)<<mean[i].imag()<<setw(26)<<errorBar[i]<<"\n";
    output_file.close();

    return 0;

}
