#include<iostream>
#include<fstream>

#include <LBFGS.h>

#include <complex.h>
#include <math.h>
#include<cmath>
#include <stdlib.h> 

#include <vector> 
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/MatrixFunctions>

using namespace Eigen;
using namespace LBFGSpp;
using namespace std;

//initialization
namespace Hsystem_ini{

   class HSystem
   {
   public: 
      int Nsite,Nspin,Ntot;

      VectorXd target;   //The amplitude of phi (Assuming it's real for simplification)
     
      HSystem(){
         Nsite=4*4;
         Nspin=3;
         
         Ntot=Nspin+Nspin;
         target=VectorXd::Zero(Nsite);

//4433u8
target(0)=   2.3042412541613433E-003;
target(1)=   4.9857982158276783E-003;
target(2)=   4.9857982165344584E-003;
target(3)=   4.9857982169589461E-003;
target(4)=   4.9857982176932346E-003;
target(5)=   1.0081847122723760E-002;
target(6)=   1.0081847125658964E-002;
target(7)=   1.0081847128341504E-002;
target(8)=   1.0081847129364723E-002;
target(9)=   1.0081847132195921E-002;
target(10)=   1.0081847134661363E-002;
target(11)=  0.49555768629773750     ;
target(12)=  0.49555768676740131     ;
target(13)= 0.49555768737821226     ;
target(14)=  0.49555768784792076     ;
target(15)=  0.93503073481453436 ;

      }
   
   };

   HSystem Hsystem;

}


// using Hsystem_ini namespace
using namespace Hsystem_ini;

//
class VMpBCS
{
private:
    int n;
public:

    void initial(VectorXd& pBCS,int i,VectorXd& mark,VectorXd& set);

    void update_mark(VectorXd& mark);

    double get_mark_value(int i,VectorXd& set,VectorXd& pBCS,VectorXd& mark);  
    
    double FNC(int Nsite, int Nspin);  
//
    double operator()(const VectorXd& v)
    {

       int total_num=int( FNC(Hsystem.Nsite,Hsystem.Nspin) );
       VectorXd pBCS=v;

       for(int i=1-1; i<=Hsystem.Nsite-1; i++){
          pBCS(i)=pBCS(i)/abs(1-pBCS(i));
          pBCS(i)=sqrt(pBCS(i));
       }

       VectorXd green= VectorXd::Zero(Hsystem.Nsite);
       VectorXd mark= VectorXd::Zero(Hsystem.Nsite-1);
       VectorXd set= VectorXd::Zero(Hsystem.Nsite-1);
       double value=0;

       for(int i=1-1; i<=Hsystem.Nsite-1; i++){

          initial(pBCS,i,mark,set);
          value=get_mark_value(i,set,pBCS,mark);
          green(i)=green(i)+value;
          for(int j=2-1; j<=total_num-1;j++){
             update_mark(mark);
             value=get_mark_value(i,set,pBCS,mark);
             green(i)=green(i)+value;
          }

       }

       green=(green/green.sum())*Hsystem.Nspin;

       double distance=0;
       for(int i=1-1;i<=Hsystem.Nsite-1;i++)distance += abs(green(i)-Hsystem.target(i)) ;
       
       return distance;
    }
};

void VMpBCS::initial(VectorXd& pBCS,int i,VectorXd& mark,VectorXd& set)
{
   int counter=0-1;
   for(int j=1-1;j<=Hsystem.Nsite-1;j++){
      if(i != j){
         counter=counter+1;
         set(counter)=pBCS(j);
      }
   }
   mark= VectorXd::Zero(Hsystem.Nsite-1);

   for(int j=1-1; j<= Hsystem.Nspin-1-1;j++){
      mark(j)=1;
   }
      
}

void VMpBCS::update_mark(VectorXd& mark)
{
   int mark_down;


   for(int j=1-1;j<=Hsystem.Nsite-2-1;j++){
      if(mark(j) == 1 && mark(j+1) == 0){
         mark(j)=0;
         mark(j+1)=1;
         mark_down=j;
         break;
      }
   }

   int left_move=0;
   for(int j=1-1;j<=mark_down-1;j++){
      if (mark(j) == 0)left_move=left_move+1;
      if (mark(j) == 1)break;
   }


   for(int j=1-1;j<=mark_down-1;j++){
      if (mark(j) == 1 && left_move != 0){
         mark(j-left_move)=1;
         mark(j)=0;
      }
   }

          
}

double VMpBCS::get_mark_value(int i,VectorXd& set,VectorXd& pBCS,VectorXd& mark)
{
   double value=pBCS(i)*pBCS(i);
   for(int j=1-1; j<=Hsystem.Nsite-1-1;j++){
      if(mark(j) == 1)value=value*set(j)*set(j);
   }
   return value;         
}

double VMpBCS::FNC(int Nsite, int Nspin)
{
   double f=1;
   for(int i=1;i<=Nsite-1;i++){
      f=f*double(i);
      if(i <= Nspin-1)f=f/double(i);
      if(i <= Nsite-Nspin)f=f/double(i);
   }
   return double(f);
}




//Numerical Hamiltonian !don't access the pubilc variables used by Hamiltonian.
class numVMpBCS
{
private:
    double step=0.0001;
public:
    VMpBCS fun;
    
    double operator()(const VectorXd& v, VectorXd& vgrad)
    {
       double distance=fun(v);
       VectorXd v_change1=v;
       VectorXd v_change2=v;
    
       for(int k=1-1;k<=Hsystem.Nsite-1;k++){
          v_change1=v;
          v_change2=v;
          v_change1( k )=v_change1( k  )+step;
          v_change2( k )=v_change2( k  )-step;
          vgrad( k  )=(fun(v_change1)-fun(v_change2))/(2*step);
       }

    return distance;
    }
};






int main()
{
    LBFGSParam<double> param;
    param.max_iterations = 3;
    param.epsilon = 0.001;
    
    LBFGSSolver<double> solver(param);
    numVMpBCS fun;
    //VMpBCS fun;

    VectorXd v = Hsystem.target;

    double fw=0;

    //test
    //VectorXd v_grad = VectorXd::Zero(Hsystem.Nsite);
    //fw=fun(v);
    //test end
    int niter = solver.minimize(fun, v, fw);

    std::cout << niter << " iterations" << std::endl;
    std::cout << "v = \n" << v << std::endl;
    std::cout << "distance from target = " << fw << std::endl;

    ofstream outfile1("VMpBCS_input.inputdat");
    outfile1<< v <<endl;
    outfile1.close();

    return 0;
}
