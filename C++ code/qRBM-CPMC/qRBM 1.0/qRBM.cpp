#include <Eigen/Core>
#include <iostream>
#include <LBFGS.h>
#include <math.h>
#include <vector> 
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
      int Nsitex,Nsitey,Nsites,Nspin_up,Nspin_dn,Ntot;
      float energy,k_energy,v_energy,sumovp;
//Hamiltonian H=-t*K+onsiteU*V
      float t,onsiteU;
// qRBM
      int NoH,powNoH; //NoH: the number of hidden {h_k} and pNoH=2^NoH

      MatrixXf Green,phi0,Hzero;

      MatrixXf *phi=NULL;  
      MatrixXf **Green_sub=NULL;   
      MatrixXf **Green_subT=NULL; 
      float **ovp=NULL; 
  
      HSystem()
      { 
         t=1;
         onsiteU=0;
      
         Nsitex=2;
         Nsitey=2; 
         Nsites=Nsitex*Nsitey;

         Nspin_up=2;
         Nspin_dn=2;
         Ntot=Nspin_up+Nspin_dn;
//the number of hidden {h_k}
         NoH=1; 

         powNoH=1; //pNoH=2^NoH
         for( int i=0; i<NoH; i++ )powNoH *= 2; 

//Def MatrixXf,vector-MatrixXf
         phi=new MatrixXf[powNoH];  

         Green_sub = new MatrixXf*[powNoH];
         for( int i=0; i<powNoH; i++ )Green_sub[i] = new MatrixXf[powNoH];

         Green_subT = new MatrixXf*[powNoH];
         for( int i=0; i<powNoH; i++ )Green_subT[i] = new MatrixXf[powNoH]; 

         ovp = new float*[powNoH];
         for( int i=0; i<powNoH; i++ )ovp[i] = new float[powNoH]; 
//Def phi0   
         phi0=MatrixXf::Zero(2*Nsites,Ntot);
   //      phi0=MatrixXf::Random(2*Nsites,Ntot);


    phi0 << 1, 1, 0, 0,
            0, 1, 0, 0,
            1, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 1, 1,
            0, 0, 0, 1,
            0, 0, 1, 0,
            0, 0, 0, 0;
         

    //phi0 << -0.19673139741673984 ,  0.84338411016005188 , 0, 0,
    //     -0.46315221476500634 ,  0.18838796659623835 , 0, 0,
    //     -0.46315221476501511 ,  0.18838796659621476 , 0, 0,
    //     -0.72957303211329505 ,  -0.46660817696758677, 0, 0,
    //     0 , 0, -0.19673139741673984, 0.84338411016005188,
    //     0 , 0,-0.46315221476500634 , 0.18838796659623835,
    //     0 , 0,-0.46315221476501511 , 0.18838796659621476,
    //     0 , 0,-0.72957303211329505 ,-0.46660817696758677;


//Def Hzero for P.B.C.
         Hzero=MatrixXf::Zero(2*Nsites,2*Nsites);
//for 2 by 2, 2up 2dn         
         Hzero(0,0)=-2; Hzero(0,1)=-0; Hzero(0,2)=-0; Hzero(0,3)=-0;
         Hzero(1,0)=-0; Hzero(1,1)=-0; Hzero(1,2)=-0; Hzero(1,3)=-0;
         Hzero(2,0)=-0; Hzero(2,1)=-0; Hzero(2,2)=-0; Hzero(2,3)=-0;
         Hzero(3,0)=-0; Hzero(3,1)=-0; Hzero(3,2)=-0; Hzero(3,3)=-0;
         Hzero.block(4,4,Nsites,Nsites)=Hzero.block(0,0,Nsites,Nsites);

      }

      float binary(int Nk, int sit)
      {
         if( (Nk >> (sit-1))%2 == 1 )
            return 1;
         else
            return -1; 
      }
   
   };

   HSystem Hsystem;

}


// using Hsystem_ini namespace
using namespace Hsystem_ini;

//
class Hamiltonian
{
private:
    int n;
public:
    
    MatrixXf get_phi(const MatrixXf& w, const int Nk, const MatrixXf& phi0);
    
    MatrixXf get_green(MatrixXf* phi);

    float get_sumovp(MatrixXf* phi);

    float get_grad(int k, int i, int j, MatrixXf** Green_sub, MatrixXf** Green_subT);

    float get_energy(const MatrixXf& Green);

//
    float operator()(const VectorXf& v, VectorXf& vgrad)
    {

    MatrixXf w=MatrixXf::Zero(Hsystem.NoH*2*Hsystem.Nsites,2*Hsystem.Nsites); 
    MatrixXf grad=MatrixXf::Zero(Hsystem.NoH*2*Hsystem.Nsites,2*Hsystem.Nsites);

// map v to w
    for(int k = 1; k <= Hsystem.NoH; k++){
       for(int i = 1; i <=Hsystem.Nsites ; i++){
          for(int j = 1; j <=Hsystem.Nsites ; j++)w(i+(k-1)*2*Hsystem.Nsites-1,j-1)=v( (k-1)*2*Hsystem.Nsites*Hsystem.Nsites+(i-1)*Hsystem.Nsites+j-1 );
       }
       for(int i = Hsystem.Nsites+1; i <=2*Hsystem.Nsites ; i++){
          for(int j = Hsystem.Nsites+1; j <=2*Hsystem.Nsites ; j++)w(i+(k-1)*2*Hsystem.Nsites-1,j-1)=v( (k-1)*2*Hsystem.Nsites*Hsystem.Nsites+(i-1)*Hsystem.Nsites+j-Hsystem.Nsites-1 );
       }
    }
  

  
    for(int Nk = 0; Nk <= Hsystem.powNoH-1; Nk++)Hsystem.phi[Nk]=get_phi(w, Nk, Hsystem.phi0);

    Hsystem.sumovp=get_sumovp(Hsystem.phi);

    Hsystem.Green=get_green(Hsystem.phi);

    //cout << w << "w" << endl;
    //cout << "Green_sub[1][1]==>:"<< endl << Hsystem.Green_sub[1][1] << endl;
    //cout << "Green_sub[0][0]==>:"<< endl << Hsystem.Green_sub[0][0] << endl;
    //cout << "Green_sub[1][0]==>:"<< endl << Hsystem.Green_sub[1][0] << endl;
    //cout << "Green_sub[0][1]==>:"<< endl << Hsystem.Green_sub[0][1] << endl;
    //cout << "Green_subT[1][1]==>:"<< endl << Hsystem.Green_subT[1][1] << endl;
    //cout << "Green_subT[0][0]==>:"<< endl << Hsystem.Green_subT[0][0] << endl;
    //cout << "Green_subT[1][0]==>:"<< endl << Hsystem.Green_subT[1][0] << endl;
    //cout << "Green_subT[0][1]==>:"<< endl << Hsystem.Green_subT[0][1] << endl;
    //cout << Hsystem.Green << endl;
    //cout << Hsystem.Hzero << endl;
    
    Hsystem.energy=get_energy(Hsystem.Green);

    

    for(int k = 1; k <= Hsystem.NoH; k++){
       for(int i = 1; i <=Hsystem.Nsites ; i++){
          for(int j = 1; j <=Hsystem.Nsites ; j++)grad(i+(k-1)*2*Hsystem.Nsites-1,j-1)=get_grad(k,i,j,Hsystem.Green_sub,Hsystem.Green_subT);
       }
       for(int i = Hsystem.Nsites+1; i <=2*Hsystem.Nsites ; i++){
          for(int j = Hsystem.Nsites+1; j <=2*Hsystem.Nsites ; j++)grad(i+(k-1)*2*Hsystem.Nsites-1,j-1)=get_grad(k,i,j,Hsystem.Green_sub,Hsystem.Green_subT);
       }
    }


// map grad to vgrad
    for(int k = 1; k <= Hsystem.NoH; k++){
       for(int i = 1; i <=Hsystem.Nsites ; i++){
          for(int j = 1; j <=Hsystem.Nsites ; j++)vgrad( (k-1)*2*Hsystem.Nsites*Hsystem.Nsites+(i-1)*Hsystem.Nsites+j -1)=grad(i+(k-1)*2*Hsystem.Nsites-1,j-1);
       }
       for(int i = Hsystem.Nsites+1; i <=2*Hsystem.Nsites ; i++){
          for(int j = Hsystem.Nsites+1; j <=2*Hsystem.Nsites ; j++)vgrad( (k-1)*2*Hsystem.Nsites*Hsystem.Nsites+(i-1)*Hsystem.Nsites+j-Hsystem.Nsites -1)=grad(i+(k-1)*2*Hsystem.Nsites-1,j-1);
       }
    }
    return Hsystem.energy;
    }
};


MatrixXf Hamiltonian::get_phi(const MatrixXf& w, const int Nk, const MatrixXf& phi0)
{
    MatrixXf phi=phi0;

    for(int k = 1; k <= Hsystem.NoH; k++){ 
       //diagonalize w(k,:,:) => wd(:,:),wv(:,:);
       //EigenSolver<MatrixXf> es( w.block((k-1)*2*Hsystem.Nsites+0,0,2*Hsystem.Nsites,2*Hsystem.Nsites)*Hsystem.binary(Nk,k) ); 
       //MatrixXf wd = es.pseudoEigenvalueMatrix(); 
       //for(int l = 0; l <= 2*Hsystem.Nsites-1; l++)wd(l,l)=exp(wd(l,l));
       //MatrixXf wv = es.pseudoEigenvectors();
       //MatrixXf ew=wv.transpose()*wd*wv;

       MatrixXf es=w.block((k-1)*2*Hsystem.Nsites+0,0,2*Hsystem.Nsites,2*Hsystem.Nsites)*Hsystem.binary(Nk,k);
       MatrixXf ew=es.exp();

       phi=ew*phi;
       }
    return phi;
}


MatrixXf Hamiltonian::get_green(MatrixXf* phi)
{
   MatrixXf I=MatrixXf::Identity(2*Hsystem.Nsites, 2*Hsystem.Nsites);
   MatrixXf Green=MatrixXf::Zero(2*Hsystem.Nsites, 2*Hsystem.Nsites);

   for(int N_k1=0; N_k1 <= Hsystem.powNoH-1; N_k1++){
      for(int N_k2=0; N_k2 <= Hsystem.powNoH-1; N_k2++){
         MatrixXf A=phi[N_k1].adjoint()*phi[N_k2];
         Hsystem.Green_sub[N_k1][N_k2]=( phi[N_k2]*A.inverse()*phi[N_k1].adjoint() ).transpose()*A.determinant();
         Hsystem.Green_subT[N_k1][N_k2]=I*Hsystem.ovp[N_k1][N_k2]-Hsystem.Green_sub[N_k1][N_k2].transpose();
         Green += Hsystem.Green_sub[N_k1][N_k2];
      }
   }
   return Green;         
}


float Hamiltonian::get_sumovp(MatrixXf* phi)
{
   float sumovp=0;

   for(int N_k1=0; N_k1 <= Hsystem.powNoH-1; N_k1++){
      for(int N_k2=0; N_k2 <= Hsystem.powNoH-1; N_k2++){
         MatrixXf A=phi[N_k1].adjoint()*phi[N_k2];
         Hsystem.ovp[N_k1][N_k2]=A.determinant();
         sumovp += Hsystem.ovp[N_k1][N_k2];
         //cout << "sumovp:" << N_k1 << "," << N_k2 << "==>:" << A.determinant() << endl;
      }
   }
   return sumovp;            
}



float Hamiltonian::get_energy(const MatrixXf& Green)
{
//K
   Hsystem.k_energy=0;
   for(int sitei=1; sitei <= 2*Hsystem.Nsites; sitei++){
      for(int sitej=1; sitej <= 2*Hsystem.Nsites; sitej++){
         Hsystem.k_energy += Hsystem.t*Green(sitei-1,sitej-1)*Hsystem.Hzero(sitei-1,sitej-1);
      }
   }
//V
   Hsystem.v_energy = 0;

//Normalization
   Hsystem.energy=(Hsystem.k_energy + Hsystem.v_energy)/Hsystem.sumovp;

}



float Hamiltonian::get_grad(int k, int i, int j, MatrixXf** Green_sub, MatrixXf** Green_subT)
{
// spin
   int opSpin=0,Spin=0,opSpin_sit=0,opSpin_sitk=0,opSpin_sitl=0;
      if( i <= Hsystem.Nsites ){
         Spin=0;
         opSpin=1;
      }      
      else{
         Spin=1;
         opSpin=-1;
      }

// d N
   float dN=0;

   for(int N_k1=0; N_k1 <= Hsystem.powNoH-1; N_k1++){
      for(int N_k2=0; N_k2 <= Hsystem.powNoH-1; N_k2++){
         dN += Hsystem.binary(N_k1,k)*Green_sub[N_k1][N_k2](i-1,j-1) + Hsystem.binary(N_k2,k)*Green_sub[N_k1][N_k2](j-1,i-1);
      }
   }
// d H
   float dK=0;
   float dV=0;

   for(int N_k1=0; N_k1 <= Hsystem.powNoH-1; N_k1++){
      for(int N_k2=0; N_k2 <= Hsystem.powNoH-1; N_k2++){
         for(int sitel=1+Spin*Hsystem.Nsites; sitel <= Hsystem.Nsites+Spin*Hsystem.Nsites; sitel++){
            for(int sitek=1+Spin*Hsystem.Nsites; sitek <= Hsystem.Nsites+Spin*Hsystem.Nsites; sitek++){                
               if(Hsystem.Hzero(sitek-1,sitel-1) == -2){
               opSpin_sitk=sitek+opSpin*Hsystem.Nsites;
               opSpin_sitl=sitel+opSpin*Hsystem.Nsites;
               //dK=0;         
                  dK += Hsystem.Hzero(sitek-1,sitel-1)*( Hsystem.binary(N_k1,k)*Green_sub[N_k1][N_k2](j-1,sitel-1)*Green_subT[N_k1][N_k2](i-1,sitek-1) + Hsystem.binary(N_k2,k)*Green_sub[N_k1][N_k2](sitek-1,j-1)*Green_subT[N_k1][N_k2](sitel-1,i-1) )/Hsystem.ovp[N_k1][N_k2];
                  dK += Hsystem.Hzero(sitek-1,sitel-1)*( Hsystem.binary(N_k1,k)*Green_sub[N_k1][N_k2](j-1,i-1)*Green_sub[N_k1][N_k2](sitek-1,sitel-1) + Hsystem.binary(N_k2,k)*Green_sub[N_k1][N_k2](i-1,j-1)*Green_sub[N_k1][N_k2](sitek-1,sitel-1) )/Hsystem.ovp[N_k1][N_k2];
                  dK += Hsystem.Hzero(sitek-1,sitel-1)*( Hsystem.binary(N_k1,k)*Green_sub[N_k1][N_k2](j-1,i-1)*Green_sub[N_k1][N_k2](opSpin_sitk-1,opSpin_sitl-1) + Hsystem.binary(N_k2,k)*Green_sub[N_k1][N_k2](i-1,j-1)*Green_sub[N_k1][N_k2](opSpin_sitk-1,opSpin_sitl-1) )/Hsystem.ovp[N_k1][N_k2];
                  //cout << "k,i,j,N_k1,N_k2,sitek,sitel "<< k << i << j << N_k1 << N_k2 << sitek << sitel << "==>" << dK <<endl;
                  //cout << "Hsystem.ovp[N_k1][N_k2] " << Hsystem.ovp[N_k1][N_k2]  << endl;
               } 
            }
         }

         //dV has some bug, let onsiteU = 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!   
         for(int sitek=1+Spin*Hsystem.Nsites; sitek <= Hsystem.Nsites+Spin*Hsystem.Nsites; sitek++){
            opSpin_sit=sitek+opSpin*Hsystem.Nsites;
            dV += Hsystem.onsiteU*( Hsystem.binary(N_k1,k)*Green_sub[N_k1][N_k2](opSpin_sit-1,opSpin_sit-1)*Green_sub[N_k1][N_k2](j-1,sitek-1)*Green_subT[N_k1][N_k2](i-1,sitek-1) + Hsystem.binary(N_k2,k)*Green_sub[N_k1][N_k2](opSpin_sit-1,opSpin_sit-1)*Green_sub[N_k1][N_k2](sitek-1,j-1)*Green_subT[N_k1][N_k2](sitek-1,i-1) );
         }
         dV +=  Hsystem.binary(N_k1,k)*Green_sub[N_k1][N_k2](j-1,i-1)*Hsystem.v_energy + Hsystem.binary(N_k2,k)*Green_sub[N_k1][N_k2](i-1,j-1)*Hsystem.v_energy;
      }
   }

   float dH = dK + dV;
   float dE = dH/Hsystem.sumovp - (dN/Hsystem.sumovp) * Hsystem.energy;
   //cout << dH << "  k .vs n  " <<dN * Hsystem.energy<< "  Hsystem.sumovp " << Hsystem.sumovp << endl;

   return dE;
}




int main()
{
    LBFGSParam<float> param;
    param.max_iterations = 1;
    //param.epsilon = 1;
    
    LBFGSSolver<float> solver(param);
    Hamiltonian fun;

    VectorXf v = VectorXf::Zero(Hsystem.NoH*2*Hsystem.Nsites*Hsystem.Nsites);
    VectorXf g_test = VectorXf::Zero(Hsystem.NoH*2*Hsystem.Nsites*Hsystem.Nsites);

    float fw=0;
    v(1)=1;
    v(4)=1;
    float fw1=fun(v,g_test);
    std::cout << "f(w1) = " << fw1 << std::endl;
    std::cout << "g_test = " << g_test << std::endl;

    v(1)=1.001;
    v(4)=1;
    float fw2=fun(v,g_test);
    std::cout << "f(w2) = " << fw2 << "<==:" << (fw2-fw1)/0.001 << std::endl;

    //v(0)=0.999;
    //float fw3=fun(v,g_test);
    //std::cout << "f(w3) = " << fw3 << "<==:" << (fw3-fw1)/-0.001 << std::endl;

    v(1)=1;
    v(4)=1.00001;
    float fw4=fun(v,g_test);
    std::cout << "f(w4) = " << fw4 << "<==:" << (fw4-fw1)/0.00001 << std::endl;

    v(1)=1.001;
    v(4)=1.001;
    float fw5=fun(v,g_test);
    std::cout << "f(w5) = " << fw5 << "<==:" << (fw5-fw1)/0.001 << std::endl;

    //int niter = solver.minimize(fun, v, fw);

    //std::cout << niter << " iterations" << std::endl;
    std::cout << "x = \n" << v.transpose() << std::endl;
    std::cout << "f(w) = " << fw << std::endl;

    return 0;
}
