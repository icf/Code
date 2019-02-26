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
         Nspin=7;
         
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

//4433u8 t0.01t0.02
target(0)=   2.3047360190699140E-003;
target(1)=   4.9687315076646678E-003;
target(2)=   4.9798367561196847E-003;
target(3)=   4.9973637749746768E-003;
target(4)=   5.0037653317879125E-003;
target(5)=   9.9778572858311744E-003;
target(6)=   1.0048959477296589E-002;
target(7)=   1.0082967971097020E-002;
target(8)=   1.0104375649555835E-002;
target(9)=   1.0120340369560648E-002;
target(10)=  1.0191952546723332E-002;
target(11)=  0.49530287841117793;     
target(12)=  0.49542307337813962;     
target(13)=  0.49567953412719967;     
target(14)=  0.49581733122279181;     
target(15)=  0.93499629617096369;

//4433u8 t0.01t0.02 pinning 0.1
target(0)=   2.3015934782632526E-003;
target(1)=   4.9593439516630016E-003;
target(2)=   4.9716478998392205E-003;
target(3)=   4.9892948426129508E-003;
target(4)=   4.9945004930035813E-003;
target(5)=    9.9199071856902361E-003;
target(6)=   9.9600103979486532E-003;
target(7)=   9.9739741341896507E-003;
target(8)=   1.0174015131442041E-002;
target(9)=   1.0178063723799209E-002;
target(10)=  1.0214687091310399E-002;
target(11)=  0.49530226356172502;    
target(12)=  0.49543869106856586;      
target(13)=  0.49569431240271261;      
target(14)=  0.49581553576199994;     
target(15)=  0.93511215887531818; 

//4477u8 t0.01t0.02
target(0)= 8.29E-02;
target(1)= 0.133400814;
target(2)= 0.133445109;
target(3)= 0.134063253;
target(4)= 0.134255672;
target(5)= 0.162702734;
target(6)= 0.163801463;
target(7)= 0.445356094;
target(8)= 0.481023445;
target(9)= 0.518402022;
target(10)= 0.551010492;
target(11)= 0.792227643;
target(12)= 0.79428524;
target(13)= 0.797388704;
target(14)= 0.79839156;
target(15)= 0.877301091;

//4455u8
target(0)=.02E-02;
target(1)=3.03E-02;
target(2)=3.05E-02;
target(3)=3.06E-02;
target(4)=3.25E-02;
target(5)=7.89E-02;
target(6)=8.06E-02;
target(7)=8.14E-02;
target(8)=8.16E-02;
target(9)=8.24E-02;
target(10)=8.42E-02;
target(11)=0.857759164;
target(12)=0.85888221;
target(13)=0.861025799;
target(14)=0.862062284;
target(15)=0.916940753;

//4477u8 t0.01t0.02
target(0)= 8.29E-02;
target(1)= 0.133400814;
target(2)= 0.133445109;
target(3)= 0.134063253;
target(4)= 0.134255672;
target(5)= 0.162702734;
target(6)= 0.163801463;
target(7)= 0.445356094;
target(8)= 0.481023445;
target(9)= 0.518402022;
target(10)= 0.551010492;
target(11)= 0.792227643;
target(12)= 0.79428524;
target(13)= 0.797388704;
target(14)= 0.79839156;
target(15)= 0.877301091;

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
       pBCS=double(Hsystem.Nspin)*pBCS/pBCS.sum();

       for(int i=1-1; i<=Hsystem.Nsite-1; i++){
          if(pBCS(i) < 0.000001){
             pBCS(i)=0.000001;
          }
          if(pBCS(i) > 0.999999999){
             pBCS(i)=0.999999999;
          }
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

       std::cout << "distance from target = " << distance << std::endl;

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
    param.max_iterations = 10;
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

    ofstream outfile1("VMpBCS_input.inputdat");
    
    v=double(Hsystem.Nspin)*v/v.sum();
    for(int i=1-1;i<=Hsystem.Nsite-1;i++){
          if(v(i) < 0.000001){
             v(i)=0.000001;
          }
          if(v(i) > 0.999999999){
             v(i)=0.999999999;
          }
        outfile1<< v(i) <<endl;
    }
    outfile1.close();

    std::cout << niter << " iterations" << std::endl;
    std::cout << "v = \n" << v << std::endl;
    std::cout << "distance from target = " << fw << std::endl;

    return 0;
}
