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
      int Dimen,Nsite,Nspin,Ntot,BC,FE_cho,Phase_lock,Dec_method,Eig_cho;
      complex<double> onsiteU,kin;
      Vector2i Nl;
      Vector2cd kbound;

      VectorXd ini_phase;

      MatrixXcd FE,Hzero,cicj_global;
      MatrixXcd eigen_values,eigen_vectors_up,eigen_vectors_dn,cicj_global_up,cicj_global_dn; // the green function (density matrix) input
    
    
      //int coor(int i, int j);

      void get_eiegns(MatrixXcd& cicj_global,MatrixXcd& eigen_values,MatrixXcd& eigen_vectors);
         
      HSystem(){
         Dimen=2;
         Nl(0)=4;Nl(1)=4;Nsite=4*4;
         kbound(0)=0.01;
         kbound(1)=0.02;
         Nspin=7;
         onsiteU=8;
         kin=-1.0;
         BC=2; //1 for PBC, 2 for TBC, 3 for open
         FE_cho=1; //1 for FE, 2 for cc_decomposition (Here we assuming FE_up=FE_dn which is ok for spin symmetry), 3 for 1_define matrix
         Eig_cho=1; //1 for eigen_values from DM, 2 for VM_pBCS (read from file)
         Phase_lock=1; // 1 for phase(Nsite-1)=0
         Dec_method=2; // 1 for DET; 2 for Analytical. 
  
         Ntot=Nspin+Nspin;

         ini_phase=VectorXd::Zero(Nsite);
         //for(int i=1-1; i<=Nsite-1; i++)ini_phase(i)=2.0*rand()/double(RAND_MAX)-1.0;
         if(Phase_lock == 1)ini_phase(Nsite-1)=0.0;


         FE=MatrixXcd::Zero(Nsite,Nspin);
         Hzero=MatrixXcd::Zero(Nsite,Nsite);   

         cicj_global=MatrixXcd::Zero(2*Nsite,2*Nsite);
         eigen_values=MatrixXcd::Zero(Nsite,Nsite);
         eigen_vectors_up=MatrixXcd::Zero(Nsite,Nsite);
         eigen_vectors_dn=MatrixXcd::Zero(Nsite,Nsite);

         cicj_global=read_cicj(); //read cicj_global from inputdat

         cicj_global_up=cicj_global.block(0,0,Nsite,Nsite);
         cicj_global_dn=cicj_global.block(Nsite,Nsite,Nsite,Nsite);

         get_eiegns(cicj_global_up,eigen_values,eigen_vectors_up); //get eigens from cicj_global
         get_eiegns(cicj_global_dn,eigen_values,eigen_vectors_dn);

         //for(int i=1-1;i<=Nsite-1;i++)std::cout  << eigen_values(i,i) << std::endl;

         if(Eig_cho == 2){
            ifstream infile1("VMpBCS_input.inputdat");
            for(int i=1-1;i<=Nsite-1;i++){
                  infile1>> eigen_values(i,i);
            }   
            infile1.close();
         }

         //test
         //eigen_vectors_up=MatrixXcd::Zero(Nsite,Nsite);
         //for(int i=1-1;i<=Nsite-1;i++)eigen_vectors_up(i,i)=1.0;

         //std::cout << "difference" << eigen_vectors_up*eigen_vectors_up.adjoint() << std::endl;

         //eigen_vectors_dn=MatrixXcd::Zero(Nsite,Nsite);
         //for(int i=1-1;i<=Nsite-1;i++)eigen_vectors_dn(i,i)=1.0;
         //end test
         Hzero=initial_lattice(BC, FE); //get the H0 and FE wave function   !!!we forced spin symmetry for Hzero!!!

 

      }

      MatrixXcd read_cicj(void)
      {
         MatrixXcd cicj_global=MatrixXcd::Zero(2*Nsite,2*Nsite);
         ifstream infile1("hubb_gdmc__obdm.dat");
            for(int i=1-1;i<=2*Nsite-1;i++){
               for(int j=1-1;j<=2*Nsite-1;j++){
                  infile1>> cicj_global(i,j);
               }
            }   
         infile1.close();
         //std::cout << "difference" << cicj_global << std::endl;
         
         return cicj_global;
      }
    

      int coor(int i, int j)   //0 is the original point
      {
         int c=-1;
         if(j == 0)c=i%Nl(0);
         if(j == 1)c=i/Nl(0);
         return c;
      }


      MatrixXcd initial_lattice(int BC, MatrixXcd& FE)
      {
          MatrixXcd Hzero,values,vectors;

          values=MatrixXcd::Zero(Nsite,Nsite);
          vectors=MatrixXcd::Zero(Nsite,Nsite);
          Hzero=MatrixXcd::Zero(Nsite,Nsite);

          if(BC == 1)Hzero=sethopPBC();
          if(BC == 2)Hzero=sethopTBC();
          if(BC == 3)Hzero=sethopOBC();


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
                          FE(k,i)=vectors(k,temp);   //check the values to make sure u get correct FE 
                          //std::cout << "dd" << vectors(k,temp) << std::endl;
                      }
                   }
                }
             }
          }  
          } // FE_cho=1 (FE=Free electrons) end

          if(FE_cho == 2){
             for(int i=1-1;i<=Nspin-1;i++){
                for(int j=1-1;j<=Nsite-1;j++)FE(j,Nspin-1-i)=eigen_vectors_up(j,Nsite-i-1);   
                //for(int j=1-1;j<=Nsite-1;j++)std::cout << "dd" << eigen_vectors_up(j,Nsite-i-1) << std::endl;
             }  
          }
      

          //test
          if(FE_cho == 3){
          FE==MatrixXcd::Zero(Nsite,Nspin);
          for(int i=1-1;i<=Nspin-1;i++)FE(i,i)=1;
          }
          //end test          

          //for(int i=1-1;i<=Nsite-1;i++)std::cout << "real_values" << real_values[i] << std::endl;
          //for(int i=1-1;i<=Nspin-1;i++)std::cout << "temp set" << temp_set[i] << std::endl;
          return Hzero;

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



   };


   void HSystem::get_eiegns(MatrixXcd& cicj_global,MatrixXcd& eigen_values,MatrixXcd& eigen_vectors)
   {
       ComplexEigenSolver<MatrixXcd> es( cicj_global ); 
       VectorXcd eigen_values_vector=es.eigenvalues();
       
       eigen_values=MatrixXcd::Zero(Nsite,Nsite);
       for(int i=1-1;i<=Nsite-1;i++)eigen_values(i,i) =  eigen_values_vector(i);
       eigen_vectors = es.eigenvectors();
       //std::cout << "dd" << eigen_values_vector << std::endl; 
       //std::cout << "dddddddddddddddddddddddddddddddddddd" << std::endl;
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

    double get_mixed_energy(MatrixXcd& pBCS, MatrixXcd& FE);

//
    double operator()(const VectorXd& v) //v is the phase of pBCS eiegn_values
    {
       MatrixXcd phase_matrix=MatrixXcd::Zero(Hsystem.Nsite,Hsystem.Nsite);
       for(int i=1-1;i<=Hsystem.Nsite-1;i++)phase_matrix(i,i)=v(i);
       MatrixXcd pBCS=get_pBCS(phase_matrix, Hsystem.eigen_values, Hsystem.eigen_vectors_up, Hsystem.eigen_vectors_dn);
 
       double energy=get_mixed_energy(pBCS,Hsystem.FE);
   
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
    
       //for(int i=1-1;i<=Hsystem.Nsite-1;i++)std::cout << "dd" << phase_dia(i,i) << std::endl;
       //for(int i=1-1;i<=Hsystem.Nsite-1;i++)std::cout << "gg" << sqrt(dia_matrix(i,i)/(1.0-dia_matrix(i,i))) << std::endl;
       //std::cout << "ddddddddddddddddddddddddddddddddddddddd" << std::endl;

       MatrixXcd pBCS=phi_up*phase_dia*phi_dn.transpose();
       
       return pBCS;
}

double VMpBCS_phase::get_mixed_energy(MatrixXcd& pBCS, MatrixXcd& FE)
{
       MatrixXcd FE_up=FE; MatrixXcd FE_dn=FE; //with spin symmetry, FE=FE_up=FE_dn
 
       MatrixXcd ovp_matrix=FE_up.transpose()*pBCS.conjugate()*FE_dn;   //<pBCS|FE>
 
       MatrixXcd Green_up=pBCS.conjugate()*FE_dn*ovp_matrix.inverse()*FE_up.transpose();
 
       MatrixXcd Green_dn=(pBCS.conjugate()).transpose()*FE_up*(ovp_matrix.inverse()).transpose()*FE_dn.transpose();

       MatrixXcd Green_pp=pBCS.conjugate()-pBCS.conjugate()*FE_dn*ovp_matrix.inverse()*FE_up.transpose()*pBCS.conjugate();     //<non projected BCS|c+_{j,up} c+_{l,dn}|W>

       MatrixXcd Green_nn=-1.0*FE_up*(ovp_matrix.inverse()).transpose()*FE_dn.transpose();   

//K
   complex<double> k_energy=0;
   for(int sitei=1-1; sitei <= Hsystem.Nsite-1; sitei++){
      for(int sitej=1-1; sitej <= Hsystem.Nsite-1; sitej++){
         k_energy += Green_up(sitei,sitej)*Hsystem.Hzero(sitei,sitej);
      }
   }

   for(int sitei=1-1; sitei <= Hsystem.Nsite-1; sitei++){
      for(int sitej=1-1; sitej <= Hsystem.Nsite-1; sitej++){
         k_energy += Green_dn(sitei,sitej)*Hsystem.Hzero(sitei,sitej);
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
  
   //calculate v_energy2 with another method
   //complex<double> v_energy2_test=0;
   //for(int i=1-1;i<=Hsystem.Nspin-1;i++){
   //   for(int j=Hsystem.Nspin+1-1;j<=Hsystem.Nsite-1;j++){
   //   v_energy2_test += conj( phase_dia(j,j)/phase_dia(i,i) );
   //   }
   //}

   //for(int i=1-1;i<=Hsystem.Nsite-1;i++)std::cout << "phase_dia" <<phase_dia(i,i)<< std::endl; 

   //MatrixXcd Green_nn_lambda=MatrixXcd::Zero(Hsystem.Nsite,Hsystem.Nsite);
   //for(int i=Hsystem.Nsite-Hsystem.Nspin+1-1;i<=Hsystem.Nsite-1;i++)Green_nn_lambda(i,i)=-1.0*( 1.0/conj(phase_dia(i,i)) );
   //std::cout << "gggggggggggggggggg" << Green_nn-Hsystem.eigen_vectors_up*Green_nn_lambda*Hsystem.eigen_vectors_dn.transpose() << std::endl;  

   //MatrixXcd Green_pp_lambda=MatrixXcd::Zero(Hsystem.Nsite,Hsystem.Nsite);
   //for(int i=1-1;i<=Hsystem.Nsite-Hsystem.Nspin-1;i++)Green_pp_lambda(i,i)=phase_dia(i,i);
   //std::cout << "gggggggggggggggggg" << Green_pp-Hsystem.eigen_vectors_up.conjugate()*Green_pp_lambda.conjugate()*Hsystem.eigen_vectors_dn.adjoint() << std::endl; 

   complex<double> sum_T=0.0;
   for(int i=1-1;i<=Hsystem.Nsite-1;i++){
      for(int j=Hsystem.Nsite-Hsystem.Nspin+1-1;j<=Hsystem.Nsite-1;j++){
         for(int k=1-1;k<=Hsystem.Nsite-Hsystem.Nspin-1;k++)sum_T += Hsystem.onsiteU*conj(phase_dia(k,k)/phase_dia(j,j))*Hsystem.eigen_vectors_up(i,j)*Hsystem.eigen_vectors_dn(i,j)*conj(Hsystem.eigen_vectors_up(i,k))*conj(Hsystem.eigen_vectors_dn(i,k));
      }
   } 
   
   std::cout << "gggggggggggggggggg: " <<energy<< "      "<< k_energy << "    "<< v_energy1 <<"     "<< v_energy2<< "  sum_T:  "<< sum_T <<"  ovp:  " << ovp_matrix.determinant() << std::endl; 
////////////////////////////////////////////
      sum_T=0.0;
      for(int j=Hsystem.Nsite-Hsystem.Nspin+1-1;j<=Hsystem.Nsite-1;j++){
         for(int k=1-1;k<=Hsystem.Nsite-Hsystem.Nspin-1;k++){ 
            complex<double> sum_temp=0.0;
            for(int i=1-1;i<=Hsystem.Nsite-1;i++)sum_temp += (Hsystem.eigen_vectors_up(i,j)*Hsystem.eigen_vectors_dn(i,j)*conj(Hsystem.eigen_vectors_up(i,k))*conj(Hsystem.eigen_vectors_dn(i,k)));
            sum_T += Hsystem.onsiteU* abs(phase_dia(k,k)/phase_dia(j,j)) *abs(sum_temp);
         }
      }
      std::cout << "real limit: " << k_energy + v_energy1 - sum_T << std::endl;   


   //diag pp and nn

   //ComplexEigenSolver<MatrixXcd> es_pp( Green_pp ); 
   //VectorXcd pp_eigen_values_vector=es_pp.eigenvalues();

   //ComplexEigenSolver<MatrixXcd> es_nn( Green_nn ); 
   //VectorXcd nn_eigen_values_vector=es_nn.eigenvalues();

   //for(int sitei=1-1; sitei <= Hsystem.Nsite-1; sitei++)std::cout << "pp" <<pp_eigen_values_vector(sitei)<< "   nn   "<< nn_eigen_values_vector(sitei) << std::endl; 

   //for(int sitei=1-1; sitei <= Hsystem.Nsite-1; sitei++)std::cout << "es_nn.eigenvectors()" << es_nn.eigenvectors().row(sitei) << std::endl; 

   //std::cout << "pp" << Hsystem.eigen_vectors_up.transpose()*Green_pp*Hsystem.eigen_vectors_up << std::endl; 

   //std::cout << "nn" << Hsystem.eigen_vectors_up.adjoint()*Green_nn*Hsystem.eigen_vectors_up.conjugate()<< std::endl;
   
   //std::cout << "A-1" <<ovp_matrix.inverse() << std::endl;

   //for(int sitei=1-1; sitei <= Hsystem.Nsite-1; sitei++)std::cout << "pBCS" << << std::endl; 
   
   return energy;

}


//Numerical Hamiltonian !don't access the pubilc variables used by Hamiltonian.
class numVMpBCS_phase
{
private:
    double step=0.001;
public:
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






int main()
{
    LBFGSParam<double> param;
    param.max_iterations = 200;
    param.epsilon = 0.001;
    
    LBFGSSolver<double> solver(param);

    numVMpBCS_phase fun;

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
    std::cout << "mixed energy with pBCS and FE = " << fw << std::endl;

    return 0;
}
