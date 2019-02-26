#define _USE_MATH_DEFINES
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

#include <algorithm>

using namespace Eigen;
using namespace LBFGSpp;
using namespace std;

//initialization
namespace Hsystem_ini{

   class HSystem
   {
   public: 
      int Dimen,Nsite,Nspin,Ntot,BC,FE_cho,Phase_lock,Dec_method,Eig_cho,ptype,Min_cho; 
      double Tot_mark_value;

      complex<double> onsiteU,kin,pinning,t_prime;
      complex<double> real_limit;
      Vector2i Nl;
      Vector2cd kbound;

      VectorXd ini_phase;

      MatrixXcd FE,FE_up,FE_dn,Hzero,Hzero_up,Hzero_dn,cicj_global,C_kl;
      MatrixXcd eigen_values,eigen_values_up,eigen_values_dn,eigen_vectors_up,eigen_vectors_dn,cicj_global_up,cicj_global_dn; // the green function (density matrix) input

      void get_eiegns(MatrixXcd& cicj_global,MatrixXcd& eigen_values,MatrixXcd& eigen_vectors);
         
      HSystem(){

         ifstream infile1("param_vmpbcs_phase");
            infile1>>Dimen;
            infile1>>Nl(0);
            infile1>>Nl(1);
            infile1>>kbound(0);
            infile1>>kbound(1);
            infile1>>Nspin;
            infile1>>onsiteU;
            infile1>>pinning;    //0.25
            infile1>>ptype;      //pinningtype: 1.+-|-+, 3.+-|+- 
            infile1>>t_prime;    //t_prime is only for closeopen B.C. 
            infile1>>kin;
            infile1>>BC;         //1 for PBC, 2 for TBC, 3 for open, 4 for closeopen without Twist
            infile1>>FE_cho;     //1 for FE, 2 for cc_decomposition (Here we assuming FE_up=FE_dn which is ok for spin symmetry), 3 for 1_define matrix
            infile1>>Eig_cho;    //1 for eigen_values from DM, 2 for VM_pBCS (read from file)
            infile1>>Phase_lock; //1 for phase(Nsite-1)=0
            infile1>>Dec_method; //1 for DET; 2 for Analytical. 
            infile1>>Min_cho;    //1 for mixed estimator (pBCS and N.O.); 2 for Abnormal term phase cancelling without Ci; 3 for Abnormal term phase cancelling              
         infile1.close();
       
         real_limit=999;
         Nsite=Nl(0)*Nl(1);
         Ntot=Nspin+Nspin;

         ini_phase=VectorXd::Zero(Nsite);
         for(int i=1-1; i<=Nsite-1; i++)ini_phase(i)=2.0*rand()/double(RAND_MAX)-1.0;
         if(Phase_lock == 1)ini_phase(Nsite-1)=0.0;


         FE=MatrixXcd::Zero(Nsite,Nspin);FE_up=MatrixXcd::Zero(Nsite,Nspin);FE_dn=MatrixXcd::Zero(Nsite,Nspin);
         Hzero=MatrixXcd::Zero(Nsite,Nsite);   
         Hzero_up=MatrixXcd::Zero(Nsite,Nsite);
         Hzero_dn=MatrixXcd::Zero(Nsite,Nsite);

         cicj_global=MatrixXcd::Zero(2*Nsite,2*Nsite);
         eigen_values=MatrixXcd::Zero(Nsite,Nsite);
         eigen_values_up=MatrixXcd::Zero(Nsite,Nsite);
         eigen_values_dn=MatrixXcd::Zero(Nsite,Nsite);
         eigen_vectors_up=MatrixXcd::Zero(Nsite,Nsite);
         eigen_vectors_dn=MatrixXcd::Zero(Nsite,Nsite);

         if(Min_cho == 2)MatrixXcd C_kl=MatrixXcd::Zero(Nsite,Nsite);

      }

      MatrixXcd read_cicj(void)
      {
         MatrixXcd cicj_global=MatrixXcd::Zero(2*Nsite,2*Nsite);

         ifstream infile1("hubb_gdmc__obdm.dat");
         if(!infile1){
            std::cout <<"read from cicj_sc_global_input" << std::endl;
            ifstream infile2("cicj_sc_global_input.inputdat");
            for(int i=1-1;i<=2*Nsite-1;i++){
               for(int j=1-1;j<=2*Nsite-1;j++){
                  infile2>> cicj_global(j,i);       
               }
            }   
            infile2.close();  
         }
         else{    
            std::cout <<"read from hubb_gdmc__obdm.dat" << std::endl;
            for(int i=1-1;i<=2*Nsite-1;i++){
               for(int j=1-1;j<=2*Nsite-1;j++){
                  infile1>> cicj_global(j,i);       //watch out!!!! there may be some problems in output of DM(i,j) or DM(j,i)
               }
            }
         }   
         infile1.close();

         return cicj_global;
      }
    

      int coor(int i, int j)   //0 is the original point
      {
         int c=-1;
         if(j == 0)c=i%Nl(0);
         if(j == 1)c=i/Nl(0);
         return c;
      }


      void initial_lattice(int BC, MatrixXcd& FE_up, MatrixXcd& FE_dn)
      {
          MatrixXcd Hzero,values,vectors;

          values=MatrixXcd::Zero(Nsite,Nsite);
          vectors=MatrixXcd::Zero(Nsite,Nsite);
          Hzero=MatrixXcd::Zero(Nsite,Nsite);
          Hzero_up=MatrixXcd::Zero(Nsite,Nsite);
          Hzero_dn=MatrixXcd::Zero(Nsite,Nsite);

          if(BC == 1)Hzero=sethopPBC();
          if(BC == 2)Hzero=sethopTBC();
          if(BC == 3)Hzero=sethopOBC();
          if(BC == 4)Hzero=sethopcloseopenBC();


          if(FE_cho == 1){
          get_eiegns(Hzero,values,vectors);

          double real_values[Nsite-1];
          for(int i=1-1;i<=Nsite-1;i++)real_values[i]=real(values(i,i));

          sort(real_values,real_values+Nsite);

          int temp=-1;
          int temp_set[Nspin-1];
          for(int j=1-1;j<=Nspin-1;j++)temp_set[j]=-1;

          for(int i=1-1;i<=Nspin-1;i++){

             temp=-1;

             for(int j=1-1;j<=Nsite-1;j++){
                if( real_values[i] == real(values(j,j)) ){
                   int flag=0;
                   temp=j;
                   for(int k=1-1;k<=Nspin-1;k++){if(temp_set[k] == temp)flag=1;}
                   if(flag == 0){
                      temp_set[i]=temp;
                      for(int k=1-1;k<=Nsite-1;k++){
                          FE_up(k,i)=vectors(k,temp);   //check the values to make sure u get correct FE
                          FE_dn(k,i)=vectors(k,temp); 
                      }
                   }
                }
             }
          }  
          } // FE_cho=1 (FE=Free electrons) end

          if(FE_cho == 2){
             for(int i=1-1;i<=Nspin-1;i++){
                for(int j=1-1;j<=Nsite-1;j++){
                   FE_up(j,Nspin-1-i)=eigen_vectors_up(j,Nsite-i-1);  
                   FE_dn(j,Nspin-1-i)=eigen_vectors_dn(j,Nsite-i-1);
                   //std::cout << "FE_up(j,Nspin-1-i):  "<< FE_up(j,Nspin-1-i) << std::endl; 
                } 
             }  
          }
      

          //test
          if(FE_cho == 3){
          FE==MatrixXcd::Zero(Nsite,Nspin);
          for(int i=1-1;i<=Nspin-1;i++){
             FE_up(i,i)=1;
             FE_dn(i,i)=1;
             }
          }
          //end test    


          //split the symmetry
          Hzero_up=Hzero;
          Hzero_dn=Hzero;
          //pinning +-|-+
          if(abs(pinning) > 0.000001){
             if(ptype == 1){
                for(int mi=1-1;mi<=Nl(0)-1;mi++){  
                  Hzero_up(mi,mi)=Hzero(mi,mi)+pinning*(pow(-1.0,(double)(mi+1) ));
                  Hzero_dn(mi,mi)=Hzero(mi,mi)-pinning*(pow(-1.0,(double)(mi+1) ));
                }
                for(int mi=1-1;mi<=Nl(0)-1;mi++){
                  Hzero_up(mi+(Nsite-Nl(0)),mi+(Nsite-Nl(0)))=Hzero(mi+(Nsite-Nl(0)),mi+(Nsite-Nl(0)))+pinning*( pow(-1.0,(double)(1+(mi+1)) ));
                  Hzero_dn(mi+(Nsite-Nl(0)),mi+(Nsite-Nl(0)))=Hzero(mi+(Nsite-Nl(0)),mi+(Nsite-Nl(0)))-pinning*( pow(-1.0,(double)(1+(mi+1)) ));
                }
             }
             if(ptype == 3){
                for(int mi=1-1;mi<=Nl(0)-1;mi++){  
                  Hzero_up(mi,mi)=Hzero(mi,mi)+pinning;
                  Hzero_dn(mi,mi)=Hzero(mi,mi)-pinning;
                }
                for(int mi=1-1;mi<=Nl(0)-1;mi++){
                  Hzero_up(mi+(Nsite-Nl(0)),mi+(Nsite-Nl(0)))=Hzero(mi+(Nsite-Nl(0)),mi+(Nsite-Nl(0)))-pinning;
                  Hzero_dn(mi+(Nsite-Nl(0)),mi+(Nsite-Nl(0)))=Hzero(mi+(Nsite-Nl(0)),mi+(Nsite-Nl(0)))+pinning;
                }
             }

          }

      
          //show Hzero
          //for(int i=1-1;i<=Nsite-1;i++){
          //   for(int j=1-1;j<=Nsite-1;j++){
          //      std::cout << " up  "<<i<<"   "<<j<<"   "<< Hzero_up(i,j) << std::endl; 
          //   }
          //}
          //for(int i=1-1;i<=Nsite-1;i++){
          //   for(int j=1-1;j<=Nsite-1;j++){
          //      std::cout << " dn  "<<i<<"   "<<j<<"   "<< Hzero_dn(i,j) << std::endl; 
          //   }
          //}


      } 

      MatrixXcd sethopPBC(void)
      {

         Hzero=MatrixXcd::Zero(Nsite,Nsite);

         for(int i=1-1;i<=Nsite-1;i++){
            for(int j=1-1;j<=Dimen-1;j++){

               int den;int ntemp;

               if(coor(i,j) == 1-1){
                 den=1;
                 for(int k=1-1;k<=j-1;k++)den=den*Nl(k);
                 ntemp=(Nl(j)-1)*den+i;
                 Hzero(i,ntemp) += kin;
                 ntemp=den+i;
                 Hzero(i,ntemp) += kin;
               }  
               else{ 
                 if(coor(i,j) == Nl(j)-1){
                    den=1;
                    for(int k=1-1;k<=j-1;k++)den=den*Nl(k);
                    ntemp=(1-Nl(j))*den+i;
                    Hzero(i,ntemp) += kin;  
                    ntemp=i-den;
                    Hzero(i,ntemp) += kin;
                 }
                 else{
                    den=1;
                    for(int k=1-1;k<=j-1;k++)den=den*Nl(k);
                    ntemp=i+den;
                    Hzero(i,ntemp) += kin;
                    ntemp=i-den;
                    Hzero(i,ntemp) += kin;
                 }
               }
            //endif
            }
         }
         return Hzero;
       
      } //end sethopPBC


      MatrixXcd sethopTBC(void)
      {

         Hzero=MatrixXcd::Zero(Nsite,Nsite);

         for(int i=1-1;i<=Nsite-1;i++){
            for(int j=1-1;j<=Dimen-1;j++){

               int den;int ntemp;

               if(coor(i,j) == 1-1){
                 den=1;
                 for(int k=1-1;k<=j-1;k++)den=den*Nl(k);
                 ntemp=(Nl(j)-1)*den+i;
                 Hzero(i,ntemp) += kin*exp(complex <double>(0,-1) *kbound(j)*2.0*M_PI/double(Nl(j)));
                 ntemp=den+i;
                 Hzero(i,ntemp) += kin*exp(complex <double>(0,1)*kbound(j)*2.0*M_PI/double(Nl(j)));
               }  
               else{ 
                 if(coor(i,j) == Nl(j)-1){
                    den=1;
                    for(int k=1-1;k<=j-1;k++)den=den*Nl(k);
                    ntemp=(1-Nl(j))*den+i;
                    Hzero(i,ntemp) += kin*exp(complex <double>(0,1)*kbound(j)*2.0*M_PI/double(Nl(j)));  
                    ntemp=i-den;
                    Hzero(i,ntemp) += kin*exp(complex <double>(0,-1)*kbound(j)*2.0*M_PI/double(Nl(j)));
                 }
                 else{
                    den=1;
                    for(int k=1-1;k<=j-1;k++)den=den*Nl(k);
                    ntemp=i+den;
                    Hzero(i,ntemp) += kin*exp(complex <double>(0,1)*kbound(j)*2.0*M_PI/double(Nl(j)));
                    ntemp=i-den;
                    Hzero(i,ntemp) += kin*exp(complex <double>(0,-1)*kbound(j)*2.0*M_PI/double(Nl(j)));
                 }
               }
            //endif
            }
         }

         return Hzero;
       
      } //end sethoTPBC

      MatrixXcd sethopOBC(void)
      {

         Hzero=MatrixXcd::Zero(Nsite,Nsite);

         for(int i=1-1;i<=Nsite-1;i++){
            for(int j=1-1;j<=Dimen-1;j++){

               int den;int ntemp;

               if(coor(i,j) == 1-1){
                 den=1;
                 for(int k=1-1;k<=j-1;k++)den=den*Nl(k);
                 ntemp=den+i;
                 Hzero(i,ntemp) += kin;
               }  
               else{ 
                 if(coor(i,j) == Nl(j)-1){
                    den=1;
                    for(int k=1-1;k<=j-1;k++)den=den*Nl(k);
                    ntemp=i-den;
                    Hzero(i,ntemp) += kin;
                 }
                 else{
                    den=1;
                    for(int k=1-1;k<=j-1;k++)den=den*Nl(k);
                    ntemp=i+den;
                    Hzero(i,ntemp) += kin;
                    ntemp=i-den;
                    Hzero(i,ntemp) += kin;
                 }
               }
            //endif
            }
         }
         return Hzero;
       
      } //end sethopOBC

      MatrixXcd sethopcloseopenBC(void)
      {

         Hzero=MatrixXcd::Zero(Nsite,Nsite);

         for(int i=1-1;i<=Nsite-1;i++){
            for(int j=1-1;j<=Dimen-1;j++){

               int den;int ntemp;

               if(coor(i,j) == 1-1){
                 den=1;
                 for(int k=1-1;k<=j-1;k++)den=den*Nl(k);
                 if(j == 1-1){
                    ntemp=(Nl(j)-1)*den+i;
                    Hzero(i,ntemp) += kin;
                 }
                 ntemp=den+i;
                 Hzero(i,ntemp) += kin;
               }  
               else{ 
                 if(coor(i,j) == Nl(j)-1){
                    den=1;
                    for(int k=1-1;k<=j-1;k++)den=den*Nl(k);
                    if(j == 1-1){
                       ntemp=(1-Nl(j))*den+i;
                       Hzero(i,ntemp) += kin; 
                    } 
                    ntemp=i-den;
                    Hzero(i,ntemp) += kin;
                 }
                 else{
                    den=1;
                    for(int k=1-1;k<=j-1;k++)den=den*Nl(k);
                    ntemp=i+den;
                    Hzero(i,ntemp) += kin;
                    ntemp=i-den;
                    Hzero(i,ntemp) += kin;
                 }
               }
            //endif
            }
         }

//add t_prime
for(int i=1-1; i<=Nsite-1; i++){
   int x1,y1,ntemp;
   int x = coor(i, 0);
   int y = coor(i, 1);
// the hopping in the ( 1, 1 ) direction
   x1 = x + 1;
   y1 = y + 1;
   if ( x1 >= Nl(0) ){
      x1 = x1 - Nl(0);
   }
   if ( y1 >= Nl(1) ){
      y1 = y1 - Nl(1);
   }
   ntemp = x1 + (y1+1-1) * Nl(0);
   Hzero(i,ntemp)  += t_prime;
   if(y1 != (y+1)){
      Hzero(i,ntemp) = 0.0;
   }
//--------------------------------------------------------------
// the hopping in the ( 1, -1 ) direction
   x1 = x + 1;
   y1 = y - 1;
   if ( x1 >= Nl(0) ){
      x1 = x1 - Nl(0);
   }
   if ( y1 < 0 ){
      y1 = y1 + Nl(1);
   }
   ntemp = x1 + (y1+1-1) * Nl(0);
   Hzero(i,ntemp)  += t_prime;
   if(y1 != (y-1)){
      Hzero(i,ntemp) = 0.0;
   }

//--------------------------------------------------------------
// the hopping in the ( -1, 1 ) direction
   x1 = x - 1;
   y1 = y + 1;
   if ( x1 < 0 ){
      x1 = x1 + Nl(0);
   }
   if ( y1 >= Nl(1) ){
      y1 = y1 - Nl(1);
   }
   ntemp = x1 + (y1+1-1) * Nl(0);
   Hzero(i,ntemp)  += t_prime;
   if(y1 != (y+1)){
      Hzero(i,ntemp) = 0.0;
   }

//--------------------------------------------------------------
// the hopping in the ( -1, -1 ) direction
   x1 = x - 1;
   y1 = y - 1;
   if ( x1 < 0 ){
      x1 = x1 + Nl(0);
   }
   if ( y1 < 0 ){
      y1 = y1 + Nl(1);
   }
   ntemp = x1 + (y1+1-1) * Nl(0);
   Hzero(i,ntemp)  += t_prime;
   if(y1 != (y-1)){
      Hzero(i,ntemp) = 0.0;
   }

}


         return Hzero;
       
      } //end sethocloseopenPBC


 //get all C_{kl} information
   void get_C_kl(void)
   {
      int Tot_num=FNC();
      Tot_mark_value=0;
      C_kl=MatrixXcd::Zero(Nsite,Nsite);

      VectorXi mark=initial_mark();
      double mark_value=get_mark_value(mark);
      Tot_mark_value += mark_value;
      update_C_kl(C_kl, mark, mark_value);

      for(int i=2-1;i<=Tot_num-1;i++){
         update_mark(mark);
         double mark_value=get_mark_value(mark);
         Tot_mark_value += mark_value;
         update_C_kl(C_kl, mark, mark_value);
      }
      std::cout <<" mark check: "<< mark << std::endl;
    
   }  
   //
   VectorXi initial_mark(void)
   {
      VectorXi mark=VectorXi::Zero(Nsite);
      for(int i=1-1;i<=Nspin-1;i++)mark(i)=1;

      return mark; 
   } 
   //
   double get_mark_value(VectorXi& mark)
   {
      double mark_value=1;
      for(int i=1-1;i<=Nsite-1;i++){
      if(mark(i) == 1)mark_value *= abs( eigen_values(i,i)/(1.0-eigen_values(i,i)) );
      }

      //mark_value rescaling
      mark_value=mark_value/10000;
     
      return mark_value; 
   } 
   //
   void update_C_kl(MatrixXcd& C_kl, VectorXi& mark, double mark_value)
   {
      for(int k=1-1;k<=Nsite-1;k++){
         for(int l=1-1;l<=Nsite-1;l++){
            if(mark(k) == 1 && mark(l) == 0){
            C_kl(k,l) += mark_value;
            }
         }
      }
   } 
   //
   void update_mark(VectorXi& mark)
   {
      int mark_down=-1;

      for(int j=1-1; j<=Nsite-1-1; j++){
         if(mark(j) == 1 && mark(j+1) == 0){
            mark(j)=0;
            mark(j+1)=1;
            mark_down=j;
            break;
         }
      }

      int left_move=0;
      for(int j=1-1; j<=mark_down-1; j++){
         if (mark(j) == 0)left_move=left_move+1;
         if (mark(j) == 1)break;
      }

      for(int j=1-1; j<=mark_down-1; j++){
         if (mark(j) == 1 && left_move != 0){
            mark(j-left_move)=1;
            mark(j)=0;
         }
      }

   }
   //choose Nspin from Nsite
   int FNC(void)
   {
      double fr=1;
      for(int i=1;i<=Nsite;i++){
         fr=fr*double(i);
         if(i <= Nspin-1)fr=fr/double(i);
         if(i <= Nsite-Nspin)fr=fr/double(i);
      }
      int f=int(fr);
      return f;
   }



   };// end HSystem


   void HSystem::get_eiegns(MatrixXcd& cicj_global,MatrixXcd& eigen_values,MatrixXcd& eigen_vectors)
   {
       ComplexEigenSolver<MatrixXcd> es( cicj_global ); 
       VectorXcd eigen_values_vector=es.eigenvalues();
       
       eigen_values=MatrixXcd::Zero(Nsite,Nsite);
       for(int i=1-1;i<=Nsite-1;i++)eigen_values(i,i) =  eigen_values_vector(i);
       eigen_vectors = es.eigenvectors();
   }  


   HSystem Hsystem;

}


// using Hsystem_ini namespace
using namespace Hsystem_ini;

//
class VMpBCS_phase
{
private:

public:

    MatrixXcd phase_dia;

    MatrixXcd get_pBCS(MatrixXcd& phase_matrix, MatrixXcd& dia_matrix, MatrixXcd& phi_up, MatrixXcd& phi_dn);

    double get_mixed_energy(MatrixXcd& pBCS, MatrixXcd& FE_up, MatrixXcd& FE_dn);

//
    double operator()(const VectorXd& v) //v is the phase of pBCS eiegn_values
    {
       MatrixXcd phase_matrix=MatrixXcd::Zero(Hsystem.Nsite,Hsystem.Nsite);
       for(int i=1-1;i<=Hsystem.Nsite-1;i++)phase_matrix(i,i)=v(i);
       MatrixXcd pBCS=get_pBCS(phase_matrix, Hsystem.eigen_values, Hsystem.eigen_vectors_up, Hsystem.eigen_vectors_dn);
 
       double energy=get_mixed_energy(pBCS,Hsystem.FE_up,Hsystem.FE_dn);
   
       return energy;
    }
};

MatrixXcd VMpBCS_phase::get_pBCS(MatrixXcd& phase_matrix, MatrixXcd& dia_matrix, MatrixXcd& phi_up, MatrixXcd& phi_dn)
{
       phase_dia=MatrixXcd::Zero(Hsystem.Nsite,Hsystem.Nsite);
       if(Hsystem.Dec_method == 1){
          for(int i=1-1;i<=Hsystem.Nsite-1;i++){
             if(i <= Hsystem.Nsite-1-Hsystem.Nspin)phase_dia(i,i)=exp( complex <double> (0,-1.0)*phase_matrix(i,i)*M_PI ) * 0.0;
             if(i > Hsystem.Nsite-1-Hsystem.Nspin)phase_dia(i,i)=exp( complex <double> (0,-1.0)*phase_matrix(i,i)*M_PI ) * 1.0;
          }        //e^{-i*phase*pi}dia
       }
       if(Hsystem.Dec_method == 2){
          for(int i=1-1;i<=Hsystem.Nsite-1;i++)phase_dia(i,i)=exp( complex <double> (0,-1.0)*phase_matrix(i,i)*M_PI ) * sqrt(dia_matrix(i,i)/(1.0-dia_matrix(i,i)));        //e^{-i*phase*pi}dia
       }

       MatrixXcd pBCS=phi_up*phase_dia*phi_dn.transpose();
       
       return pBCS;
}

double VMpBCS_phase::get_mixed_energy(MatrixXcd& pBCS, MatrixXcd& FE_up, MatrixXcd& FE_dn)
{
       MatrixXcd ovp_matrix=FE_up.transpose()*pBCS.conjugate()*FE_dn;   //<pBCS|FE>
 
       MatrixXcd Green_up=pBCS.conjugate()*FE_dn*ovp_matrix.inverse()*FE_up.transpose();
 
       MatrixXcd Green_dn=(pBCS.conjugate()).transpose()*FE_up*(ovp_matrix.inverse()).transpose()*FE_dn.transpose();

       MatrixXcd Green_pp=pBCS.conjugate()-pBCS.conjugate()*FE_dn*ovp_matrix.inverse()*FE_up.transpose()*pBCS.conjugate();     //<non projected BCS|c+_{j,up} c+_{l,dn}|W>

       MatrixXcd Green_nn=-1.0*FE_up*(ovp_matrix.inverse()).transpose()*FE_dn.transpose();   

//K
   complex<double> k_energy=0;
   for(int sitei=1-1; sitei <= Hsystem.Nsite-1; sitei++){
      for(int sitej=1-1; sitej <= Hsystem.Nsite-1; sitej++){
         k_energy += Green_up(sitei,sitej)*Hsystem.Hzero_up(sitei,sitej);
      }
   }

   for(int sitei=1-1; sitei <= Hsystem.Nsite-1; sitei++){
      for(int sitej=1-1; sitej <= Hsystem.Nsite-1; sitej++){
         k_energy += Green_dn(sitei,sitej)*Hsystem.Hzero_dn(sitei,sitej);
      }
   }

//V
   complex<double> v_energy1=0;
   complex<double> v_energy2=0;
   for(int sitei=1-1; sitei <= Hsystem.Nsite-1; sitei++){
      v_energy1 += Hsystem.onsiteU*Green_up(sitei,sitei)*Green_dn(sitei,sitei);
      v_energy2 -= Hsystem.onsiteU*Green_pp(sitei,sitei)*Green_nn(sitei,sitei);
   }

   double energy=real(k_energy + v_energy1 + v_energy2);

   complex<double> sum_T=0.0;
   for(int i=1-1;i<=Hsystem.Nsite-1;i++){
      for(int j=Hsystem.Nsite-Hsystem.Nspin+1-1;j<=Hsystem.Nsite-1;j++){
         for(int k=1-1;k<=Hsystem.Nsite-Hsystem.Nspin-1;k++)sum_T += Hsystem.onsiteU*conj(phase_dia(k,k)/phase_dia(j,j))*Hsystem.eigen_vectors_up(i,j)*Hsystem.eigen_vectors_dn(i,j)*conj(Hsystem.eigen_vectors_up(i,k))*conj(Hsystem.eigen_vectors_dn(i,k));
      }
   } 

    //std::cout << "gggggggggggggggggg: " <<energy<< "      "<< k_energy << "    "<< v_energy1 <<"     "<< v_energy2<< "  sum_T:  "<< sum_T <<"  ovp:  " << ovp_matrix.determinant() << std::endl; 

      sum_T=0.0;
      for(int j=Hsystem.Nsite-Hsystem.Nspin+1-1;j<=Hsystem.Nsite-1;j++){
         for(int k=1-1;k<=Hsystem.Nsite-Hsystem.Nspin-1;k++){ 
            complex<double> sum_temp=0.0;
            for(int i=1-1;i<=Hsystem.Nsite-1;i++)sum_temp += (Hsystem.eigen_vectors_up(i,j)*Hsystem.eigen_vectors_dn(i,j)*conj(Hsystem.eigen_vectors_up(i,k))*conj(Hsystem.eigen_vectors_dn(i,k)));
            sum_T += Hsystem.onsiteU* abs(phase_dia(k,k)/phase_dia(j,j)) *abs(sum_temp);
         }
      }
    Hsystem.real_limit=k_energy + v_energy1 - sum_T;
    //std::cout << "real limit: " << k_energy + v_energy1 - sum_T << std::endl;  


    if(Hsystem.Min_cho == 2){
      MatrixXcd M_jk=MatrixXcd::Zero(Hsystem.Nsite,Hsystem.Nsite);
      for(int j=1-1;j<=Hsystem.Nsite-1;j++){
         for(int k=1-1;k<=Hsystem.Nsite-1;k++){ 
            complex<double> sum_temp=0.0;
            for(int i=1-1;i<=Hsystem.Nsite-1;i++)sum_temp += (Hsystem.eigen_vectors_up(i,j)*Hsystem.eigen_vectors_dn(i,j)*conj(Hsystem.eigen_vectors_up(i,k))*conj(Hsystem.eigen_vectors_dn(i,k)));
            M_jk(j,k) = sum_temp;
         }
      }
      energy=0;
      Hsystem.real_limit=0;
      for(int j=Hsystem.Nsite-Hsystem.Nspin+1-1;j<=Hsystem.Nsite-1;j++){
         for(int k=1-1;k<=Hsystem.Nsite-Hsystem.Nspin-1;k++){ 
            energy += real( conj( phase_dia(k,k)/phase_dia(j,j) )*M_jk(j,k) );
            energy += real( conj( phase_dia(j,j)/phase_dia(k,k) )*M_jk(k,j) );
            Hsystem.real_limit += -1.0*abs( conj( phase_dia(k,k)/phase_dia(j,j) )*M_jk(j,k) );
            Hsystem.real_limit += -1.0*abs( conj( phase_dia(j,j)/phase_dia(k,k) )*M_jk(k,j) );
         }
      }
      std::cout <<" energy(abnormal): "<< energy <<" real limit: " << Hsystem.real_limit << std::endl;  

      energy=0;
      Hsystem.real_limit=0;
      for(int j=1-1;j<=Hsystem.Nsite-1;j++){
         for(int k=1-1;k<=Hsystem.Nsite-1;k++){ 
            energy += real( conj( phase_dia(k,k)/phase_dia(j,j) )*M_jk(j,k) );
            Hsystem.real_limit += -1.0*abs( conj( phase_dia(k,k)/phase_dia(j,j) )*M_jk(j,k) );
         }
      }
      std::cout <<" total energy(abnormal): "<< energy <<" total real limit: " << Hsystem.real_limit << std::endl;  

    }


    if(Hsystem.Min_cho == 3){
      MatrixXcd M_jk=MatrixXcd::Zero(Hsystem.Nsite,Hsystem.Nsite);
      for(int j=1-1;j<=Hsystem.Nsite-1;j++){
         for(int k=1-1;k<=Hsystem.Nsite-1;k++){ 
            complex<double> sum_temp=0.0;
            for(int i=1-1;i<=Hsystem.Nsite-1;i++)sum_temp += (Hsystem.eigen_vectors_up(i,j)*Hsystem.eigen_vectors_dn(i,j)*conj(Hsystem.eigen_vectors_up(i,k))*conj(Hsystem.eigen_vectors_dn(i,k)));
            M_jk(j,k) = sum_temp;
         }
      }
      energy=0;
      Hsystem.real_limit=0;
      for(int j=Hsystem.Nsite-Hsystem.Nspin+1-1;j<=Hsystem.Nsite-1;j++){
         for(int k=1-1;k<=Hsystem.Nsite-Hsystem.Nspin-1;k++){ 
            energy += real( (Hsystem.C_kl(j,k)/Hsystem.Tot_mark_value)*conj( phase_dia(k,k)/phase_dia(j,j) )*M_jk(j,k) );
            energy += real( (Hsystem.C_kl(k,j)/Hsystem.Tot_mark_value)*conj( phase_dia(j,j)/phase_dia(k,k) )*M_jk(k,j) );
            Hsystem.real_limit += -1.0*abs( (Hsystem.C_kl(j,k)/Hsystem.Tot_mark_value)*conj( phase_dia(k,k)/phase_dia(j,j) )*M_jk(j,k) );
            Hsystem.real_limit += -1.0*abs( (Hsystem.C_kl(k,j)/Hsystem.Tot_mark_value)*conj( phase_dia(j,j)/phase_dia(k,k) )*M_jk(k,j) );
         }
      }
      std::cout <<" energy(abnormal): "<< energy <<" real limit: " << Hsystem.real_limit << std::endl;  

      energy=0;
      Hsystem.real_limit=0;
      for(int j=1-1;j<=Hsystem.Nsite-1;j++){
         for(int k=1-1;k<=Hsystem.Nsite-1;k++){ 
            energy += real( (Hsystem.C_kl(j,k)/Hsystem.Tot_mark_value)*conj( phase_dia(k,k)/phase_dia(j,j) )*M_jk(j,k) );
            Hsystem.real_limit += -1.0*abs( (Hsystem.C_kl(j,k)/Hsystem.Tot_mark_value)*conj( phase_dia(k,k)/phase_dia(j,j) )*M_jk(j,k) );
         }
      }
      std::cout <<" total energy(abnormal): "<< energy <<" total real limit: " << Hsystem.real_limit << std::endl;  

    }

   
   
   return energy;

}


//Numerical Hamiltonian !don't access the pubilc variables used by Hamiltonian.
class numVMpBCS_phase
{
private:
    
public:
    double step=0.001;
    VMpBCS_phase fun;
    
    double operator()(const VectorXd& v, VectorXd& vgrad)
    {
       double energy=fun(v);
       VectorXd v_change1=v;
       VectorXd v_change2=v;
    
       for(int k=1-1;k<=Hsystem.Nsite-1;k++){
          if(Hsystem.Phase_lock != 1 || k!=Hsystem.Nsite-1){
             v_change1=v;
             v_change2=v;
             v_change1( k )=v_change1( k  )+step;
             v_change2( k )=v_change2( k  )-step;
             vgrad( k  )=(fun(v_change1)-fun(v_change2))/(2.0*step);
          }
          else{
             vgrad( k )=0.0;
          }
       }

    return energy;
    }
};

extern "C"{
   void vmpbcs_phase_();
   class numVMpBCS_phase;
   class VMpBCS_phase;
   class HSystem;
}

void vmpbcs_phase_()
{

    LBFGSParam<double> param;
    param.max_iterations = 200;
    param.epsilon = 0.001;
    
    LBFGSSolver<double> solver(param);

    numVMpBCS_phase fun;

    //initial 
    Hsystem.cicj_global=Hsystem.read_cicj(); //read cicj_global from inputdat
    // make the DM hermitian
    

    MatrixXcd ad_temp=Hsystem.cicj_global.adjoint();
    Hsystem.cicj_global=(Hsystem.cicj_global+ad_temp )/2.0;  

    Hsystem.cicj_global_up=(Hsystem.cicj_global).block(0,0,Hsystem.Nsite,Hsystem.Nsite);
    Hsystem.cicj_global_dn=(Hsystem.cicj_global).block(Hsystem.Nsite,Hsystem.Nsite,Hsystem.Nsite,Hsystem.Nsite);

    Hsystem.get_eiegns(Hsystem.cicj_global_up,Hsystem.eigen_values_up,Hsystem.eigen_vectors_up); //get eigens from cicj_global
    Hsystem.get_eiegns(Hsystem.cicj_global_dn,Hsystem.eigen_values_dn,Hsystem.eigen_vectors_dn);

/////////////////////////////////////////////////////////////////////////////////////////
    ofstream outfile_va_up("VMpBCS_phase_eigenvalues_up.dat");
    for(int i=1-1;i<=Hsystem.Nsite-1;i++)outfile_va_up<< Hsystem.eigen_values_up(i,i) <<endl;
    outfile_va_up.close();
    ofstream outfile_va_dn("VMpBCS_phase_eigenvalues_dn.dat");
    for(int i=1-1;i<=Hsystem.Nsite-1;i++)outfile_va_dn<< Hsystem.eigen_values_dn(i,i) <<endl;
    outfile_va_dn.close();

/////////////////////////////////////////////////////////////////////////////////////////
    ofstream outfile_ve_up("VMpBCS_phase_eigenvectors_up.dat");
    for(int i=1-1;i<=Hsystem.Nsite-1;i++){
       for(int j=1-1;j<=Hsystem.Nsite-1;j++)outfile_ve_up<< Hsystem.eigen_vectors_up(j,i) <<endl;
    }
    outfile_ve_up.close();
    ofstream outfile_ve_dn("VMpBCS_phase_eigenvectors_dn.dat");
    for(int i=1-1;i<=Hsystem.Nsite-1;i++){
       for(int j=1-1;j<=Hsystem.Nsite-1;j++)outfile_ve_dn<< Hsystem.eigen_vectors_dn(j,i) <<endl;
    }
    outfile_ve_dn.close();

/////////////////////////////////////////////////////////////////////////////////////////

    Hsystem.eigen_values=(Hsystem.eigen_values_up+Hsystem.eigen_values_dn)/2.0;

    if(Hsystem.Eig_cho == 2){
       ifstream infile1("VMpBCS_input.inputdat");
       for(int i=1-1;i<=Hsystem.Nsite-1;i++){
             infile1>> Hsystem.eigen_values(i,i);
       }   
       infile1.close();
    }

    if(Hsystem.Min_cho == 3)Hsystem.get_C_kl();
    //std::cout <<" C_kl: "<< Hsystem.C_kl << std::endl; 

    Hsystem.initial_lattice(Hsystem.BC, Hsystem.FE_up, Hsystem.FE_dn); //get the H0 and FE wave function  

    //
    double fw=0;

    VectorXd v = Hsystem.ini_phase;
    //test
    //VectorXd v_grad = VectorXd::Zero(Hsystem.Nsite);
    //fw=fun(v,v_grad);
    //test end

    int niter = solver.minimize(fun, v, fw);

    ofstream outfile1("VMpBCS_phase_input.inputdat");
    for(int i=1-1;i<=Hsystem.Nsite-1;i++)outfile1<< v(i) <<endl;
    outfile1.close();
    std::cout << niter << " iterations" << std::endl;

    std::cout << "v = \n" << v << std::endl;
    std::cout << "mixed energy with pBCS and FE = " << fw <<" real_limit:  "<< Hsystem.real_limit<< std::endl;

}


