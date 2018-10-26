#include<iostream>
#include<fstream>
#include<iomanip>
#include<cstdlib>
int main( )
{
    std::ofstream out_file;
    out_file.open("k_point_1") ;
    srand (time(NULL));
    std::cout << time(NULL) << std::endl ;
    for( int i = 0 ; i < 8 ; i++ )
    {
        double a = ( rand( )%100 )/100.0;
        double b = ( rand( )%100 )/100.0;
        //std::cout << a << "\t" << b << std::endl ;
        out_file <<std::fixed << std::setprecision(2) << a <<"\t" << b << "\n" ;
    }
    out_file.close( ) ;
}
