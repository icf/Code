#include "../include/vmpBCS.h"

using namespace tensor_hao;
using namespace Eigen;
using namespace LBFGSpp;
using namespace std;

void GpBCS_VM::setModel(Model* model_, string initialSCPhiTFlag_)
{
    model=model_;
    initialSCPhiTFlag=initialSCPhiTFlag_;
    Nsite=model->getL();
    Ntot=model->getN();
    Nspin=Ntot/2;
    Tot_sample_num=20000;

    if( initialSCPhiTFlag == "setFromDensity_VMGpBCS"){
       N_vm=Nsite;
    }else if( initialSCPhiTFlag == "setFromGHF"){
       N_vm=2*Nsite;
    }else{
      cout<<"initialSCPhiTFlag Error1: "<<initialSCPhiTFlag<<endl;
      exit(1);
    }

} 

void GpBCS_VM::set(MatrixXcd& target_dm_, MatrixXcd& eigenvalues, MatrixXcd& eigenvectors, VectorXcd& v)
{
    target_dm=target_dm_;
    dia=MatrixXcd::Zero(2*Nsite,2*Nsite);
    savedDia=MatrixXcd::Zero(2*Nsite,2*Nsite);
    target_dia=MatrixXcd::Zero(2*Nsite,2*Nsite);
    T=MatrixXcd::Zero(2*Nsite,2*Nsite);
    T_up=MatrixXcd::Zero(Nsite,Nsite);
    T_dn=MatrixXcd::Zero(Nsite,Nsite);

    //dia
    dia=eigenvalues;
    savedDia=dia;

    updateDia(v);
    
    //T
    if( initialSCPhiTFlag == "setFromDensity_VMGpBCS"){
       for(int i=1-1;i<=2*Nsite-1;i++){
       for(int j=1-1;j<=2*Nsite-1;j++){
          T(i,j)=eigenvectors(i,j);
       }
       }
       for(int i=1-1;i<=Nsite-1;i++){
       for(int j=1-1;j<=Nsite-1;j++){
          T_up(i,j)=T(i,2*j);
          T_dn(i,j)=T(i+Nsite,2*j+1);
       }
       }
    }else if( initialSCPhiTFlag == "setFromGHF"){
       
       Ghf ghf;
       TensorHao<complex<double>, 2> T_hao=ghf.run(); 
       TensorHao<complex<double>, 2> T_flip(2*Nsite,2*Nsite);     

       for(int i=1-1;i<=2*Nsite-1;i++){
       for(int j=1-1;j<=2*Nsite-1;j++){
          T_flip(i,j)=T_hao(i,j);
       }
       }
       for(int i=1-1;i<=2*Nsite-1;i++){
       for(int j=1-1;j<=2*Nsite-1;j++){
          T(i,2*Nsite-1-j)=T_flip(i,j);
       }
       }

    }else{
      cout<<"initialSCPhiTFlag Error2: "<<initialSCPhiTFlag<<endl;
      exit(1);
    }

    //M=T*dia*T^T
    M= T*dia*T.transpose(); 
}

VectorXi GpBCS_VM::sample_initial_mark(void)
{
   srand((unsigned int)(time(NULL)));
   VectorXi mark=VectorXi::Zero(Nsite);
   for(int i=1-1;i<=Nspin-1;i++)mark(i)=1;  //Attention! make sure the last term is the largest
   return mark; 
} 

//
double GpBCS_VM::get_mark_value(VectorXi& mark)    //a fast estimation of <G|i><i|G>
{
   double mark_value=1;
   for(int i=1-1;i<=Nsite-1;i++){
      if(mark(i) == 1)mark_value *= abs(dia(2*i,2*i+1))*abs(dia(2*i,2*i+1));   //what is the value of mark_value?   we need to sample abs(dia(2*i,2*i+1))*abs(dia(2*i,2*i+1)) to get pure energy!!
   }
   return mark_value; 
} 

//
MatrixXcd GpBCS_VM::get_mark_SD(VectorXi& mark)
{
   int mark_cho=-1;     
   MatrixXcd SD_marked=MatrixXcd::Zero(2*Nsite,Ntot);
   for(int i=1-1;i<=Nsite-1;i++){
      if(mark(i) == 1){
         mark_cho++;
         SD_marked(2*i,mark_cho)=1.0;
         mark_cho++;
         SD_marked(2*i+1,mark_cho)=1.0;
      }
   }
   if(mark_cho != Ntot-1){
      cout<<"get SD_marked error: "<<mark_cho<<" standard "<<Ntot-1<<endl;
   }

   SD_marked=T*SD_marked;
   
   return SD_marked; 
} 

void GpBCS_VM::sample_update_mark(VectorXi& mark)
{
   //srand((unsigned int)(time(NULL)));
   int rand_pick=(rand()%(Nsite-Nspin-1+1))+1;
   int counter=0;
   int pick_mark_down=-1;

   for(int i=1-1;i<=Nsite-1;i++){
      if(mark(i) == 0)counter += 1;
      if(counter == rand_pick){
         pick_mark_down=i;
         break;
      }
   }

   //srand((unsigned int)(time(NULL)));
   int rand_drop=(rand()%(Nspin-1+1))+1;
   counter=0;
   int drop_mark_down=-1;

   for(int i=1-1;i<=Nsite-1;i++){
      if(mark(i) == 1)counter += 1;
      if(counter == rand_drop){
         drop_mark_down=i;
         break;
      }
   }

   VectorXi mark_try=VectorXi::Zero(Nsite);
   mark_try=mark;
   mark_try(pick_mark_down)=1;
   mark_try(drop_mark_down)=0;
   
   double value0=get_mark_value(mark);  
   double value1=get_mark_value(mark_try);
   
   if( double(rand())/double( RAND_MAX ) <= sqrt(value1/value0) )mark=mark_try;    
}

double GpBCS_VM::e_v_update(VectorXcd& v)
{
   updateDia(v);

   //M=T*dia*T^T
   M=T*dia*T.transpose();

   if( initialSCPhiTFlag == "setFromDensity_VMGpBCS"){
      //get energy from GpBCS (<GpBCS|CiCj|GpBCS>/<GpBCS|GpBCS>)
      //getEnergyFastMixed();
      //getEnergyFast();
      //getEnergyFastHFB();    //Phase won't affect <H> of HFB
      getEnergyFastICF();      //Calculate pure pBCS estimator in icf summing, ATTENTION: need to set C_jk first (before VM)
   }else if( initialSCPhiTFlag == "setFromGHF"){
      //get energy from GpBCS (<GpBCS|CiCj|GpBCS>/<GpBCS|GpBCS>)
      getEnergyFast();
   }else{
      cout<<"initialSCPhiTFlag Error3: "<<initialSCPhiTFlag<<endl;
      exit(1);
   }

   return energy;
}

void GpBCS_VM::updateDia(VectorXcd& v)  // a can be any unitary transformation
{
   if( initialSCPhiTFlag == "setFromDensity_VMGpBCS"){
      int counter=0;
      dia=MatrixXcd::Zero(2*Nsite,2*Nsite);
      for(int i=1-1;i<=Nsite-1;i++){
          dia(2*i,2*i+1)=exp(complex <double> (0,-1)*v(counter)*M_PI)*savedDia(2*i,2*i+1);
          dia(2*i+1,2*i)=-1.0*dia(2*i,2*i+1);
          counter++;
      }
      if(counter != N_vm){
         cout<<"number of vm error: "<<counter<<" standard: "<<N_vm<<endl;
      }
   }else if( initialSCPhiTFlag == "setFromGHF"){
      int counter=0;
      dia=MatrixXcd::Zero(2*Nsite,2*Nsite);
      for(int i=1-1;i<=Nsite-1;i++){
          dia(2*i,2*i+1)=exp(complex <double> (0,-1)*v(counter)*M_PI)*v(counter+1);
          dia(2*i+1,2*i)=-1.0*dia(2*i,2*i+1);
          counter++;counter++;
      }
      if(counter != N_vm){
         cout<<"number of vm error: "<<counter<<" standard: "<<N_vm<<endl;
      }
   }else{
         cout<<"initialSCPhiTFlag Error4: "<<initialSCPhiTFlag<<endl;
         exit(1);
   }

}

void GpBCS_VM::getEnergyFast()
{
   int Tot_num=Tot_sample_num;
   double mark_sum=0;

   VectorXi mark=sample_initial_mark();

   energy=get_mixed_energy(mark)*sqrt(get_mark_value(mark)); 
   mark_sum += sqrt(get_mark_value(mark));

   for(int i=2-1;i<=Tot_num-1;i++){
      sample_update_mark(mark);
      energy +=get_mixed_energy(mark)*sqrt(get_mark_value(mark));
      mark_sum += sqrt(get_mark_value(mark));
   }
   energy=energy/mark_sum;
}

void GpBCS_VM::getEnergyFastMixed()
{
   VectorXi mark=sample_initial_mark();
   energy=get_mixed_energy(mark);
}

void GpBCS_VM::getEnergyFastHFB()
{
   energy=get_HFB_energy();     //icf: not sure if <HFB||HFB> follow wick's theory
}

void GpBCS_VM::getEnergyFastICF()
{
   energy=get_pure_abEnergy_icf();     //icf: not sure if <HFB||HFB> follow wick's theory
}

void GpBCS_VM::getDensityMatrixDiatance()
{
   int Tot_num=Tot_sample_num;
   double mark_sum=0;
   double energy_save=0;
   double temp_energy=0;

   VectorXi mark=sample_initial_mark();

   dm=MatrixXcd::Zero(2*Nsite,2*Nsite);
   dm=get_mixed_dm(mark,energy_save)*sqrt(get_mark_value(mark)); 
   temp_energy=energy_save*sqrt(get_mark_value(mark));
   mark_sum += sqrt(get_mark_value(mark));

   for(int i=2-1;i<=Tot_num-1;i++){
      sample_update_mark(mark);
      dm +=get_mixed_dm(mark,energy_save)*sqrt(get_mark_value(mark));
      temp_energy +=energy_save*sqrt(get_mark_value(mark));
      mark_sum += sqrt(get_mark_value(mark));
   }
   temp_energy=temp_energy/mark_sum;

   dm=dm/mark_sum;

   double tempSum=0.0;
   for(int i=1-1;i<=2*Nsite-1;i++){
   for(int j=1-1;j<=2*Nsite-1;j++){
      tempSum += abs(dm(i,j)-target_dm(i,j))*abs(dm(i,j)-target_dm(i,j));
   }
   }

   energy=10.0*tempSum+temp_energy;
}

double GpBCS_VM::get_mixed_energy(VectorXi& mark)
{
   int L  = Nsite; int L2 = L*2; int N=Ntot;
   MatrixXcd SD=get_mark_SD(mark);

   TensorHao<complex<double>, 2> SD_hao_temp(L2,N);
   for(int i=1-1;i<=L2-1;i++){
   for(int j=1-1;j<=N-1;j++){
      SD_hao_temp(i,j)=SD(i,j);
   }
   }
   WalkerRight SD_hao;
   SD_hao.resize(L2,N);
   SD_hao.wfRef()=SD_hao_temp;
   SD_hao.logwRef()=complex<double> (0,0);

   TensorHao<complex<double>, 2> GpBCS_hao_temp(L2,L2);
   for(int i=1-1;i<=L2-1;i++){
   for(int j=1-1;j<=L2-1;j++){
      GpBCS_hao_temp(i,j)=M(i,j);
   }
   }
   BCS GpBCS_hao;
   GpBCS_hao.resize(L2);
   GpBCS_hao.wfRef()=GpBCS_hao_temp;
   GpBCS_hao.logwRef()=complex<double> (0,0);
  
   BCSSDOperation bcssdOperation;
   bcssdOperation.set(GpBCS_hao, SD_hao);
   const TensorHao< complex<double>, 2 > &greenMatrix = bcssdOperation.returnGreenMatrix();
   const TensorHao< complex<double>, 2 > &cmcmMatrix = bcssdOperation.returnCmCmMatrix();
   const TensorHao< complex<double>, 2 > &cpcpMatrix = bcssdOperation.returnCpCpMatrix();
   complex<double> Kenergy(0,0), Venergy(0,0);
   //Add K
   const TensorHao< complex<double>, 2 > &K = model->getK();
   for(int i = 0; i < L2; ++i)
   {
       for(int j = 0; j < L2; ++j)
       {
          Kenergy += K(j,i) * greenMatrix(j,i);
       }
   }

   //Add U
   const TensorHao< double, 1> &U = model->getU();
   for(int i = 0; i < L; ++i)
   {
       Venergy += U(i) * ( greenMatrix(i,i)*greenMatrix(i+L,i+L) - greenMatrix(i, i+L)*greenMatrix(i+L, i) );
       Venergy += U(i) * (-1.0* cpcpMatrix(i,i+L)*cmcmMatrix(i,i+L) );
   } 
   double energy = real( Kenergy + Venergy );
   return energy;
}

MatrixXcd GpBCS_VM::get_mixed_dm(VectorXi& mark, double &Energy)
{
   int L  = Nsite; int L2 = L*2; int N=Ntot;
   MatrixXcd SD=get_mark_SD(mark);

   TensorHao<complex<double>, 2> SD_hao_temp(L2,N);
   for(int i=1-1;i<=L2-1;i++){
   for(int j=1-1;j<=N-1;j++){
      SD_hao_temp(i,j)=SD(i,j);
   }
   }
   WalkerRight SD_hao;
   SD_hao.resize(L2,N);
   SD_hao.wfRef()=SD_hao_temp;
   SD_hao.logwRef()=complex<double> (0,0);

   TensorHao<complex<double>, 2> GpBCS_hao_temp(L2,L2);
   for(int i=1-1;i<=L2-1;i++){
   for(int j=1-1;j<=L2-1;j++){
      GpBCS_hao_temp(i,j)=M(i,j);
   }
   }
   BCS GpBCS_hao;
   GpBCS_hao.resize(L2);
   GpBCS_hao.wfRef()=GpBCS_hao_temp;
   GpBCS_hao.logwRef()=complex<double> (0,0);
  
   BCSSDOperation bcssdOperation;
   bcssdOperation.set(GpBCS_hao, SD_hao);
   const TensorHao< complex<double>, 2 > &greenMatrix = bcssdOperation.returnGreenMatrix();
   const TensorHao< complex<double>, 2 > &cmcmMatrix = bcssdOperation.returnCmCmMatrix();
   const TensorHao< complex<double>, 2 > &cpcpMatrix = bcssdOperation.returnCpCpMatrix();
   complex<double> Kenergy(0,0), Venergy(0,0);
   //Add K
   const TensorHao< complex<double>, 2 > &K = model->getK();
   for(int i = 0; i < L2; ++i)
   {
       for(int j = 0; j < L2; ++j)
       {
          Kenergy += K(j,i) * greenMatrix(j,i);
       }
   }

   //Add U
   const TensorHao< double, 1> &U = model->getU();
   for(int i = 0; i < L; ++i)
   {
       Venergy += U(i) * ( greenMatrix(i,i)*greenMatrix(i+L,i+L) - greenMatrix(i, i+L)*greenMatrix(i+L, i) );
       Venergy += U(i) * (-1.0* cpcpMatrix(i,i+L)*cmcmMatrix(i,i+L) );
   }

   Energy = real( Kenergy + Venergy );

   MatrixXcd mixed_dm=MatrixXcd::Zero(L2,L2);
   for(int i=1-1;i<=L2-1;i++){
   for(int j=1-1;j<=L2-1;j++){
      mixed_dm(i,j)=greenMatrix(i,j);
   }
   }
   
   return mixed_dm;
}

double GpBCS_VM::get_HFB_energy()
{
   int L  = Nsite; int L2 = L*2; int N=Ntot;
   //get HFB U V presentation
   TensorHao<complex<double>, 2> T_Hao(L2,L2);
   for(int i=1-1;i<=L2-1;i++){
   for(int j=1-1;j<=L2-1;j++){
      T_Hao(i,j)=T(i,j);
   }
   }
   //get lambda_v and lambda_u
   TensorHao<complex<double>, 2> lambda_v(L2,L2),lambda_u(L2,L2);lambda_v=complex <double> (0.0,0.0);lambda_u=complex <double> (0.0,0.0);
   for(size_t i=1-1;i<=L2/2-1;i++){
      complex <double> tempParam,tempLogPhase;
      tempParam=exp(log(dia(2*i,2*i+1)).real());
      tempLogPhase=log(dia(2*i,2*i+1)).imag();
      //
      lambda_v(2*i,2*i+1)=exp(complex<double>(0.0,-1.0)*tempLogPhase)*sqrt(tempParam*tempParam/(1.0+tempParam*tempParam));
      lambda_v(2*i+1,2*i)=lambda_v(2*i,2*i+1);
      //
      lambda_u(2*i,2*i)=1.0/sqrt(1.0+tempParam*tempParam);
      lambda_u(2*i+1,2*i+1)=-1.0*lambda_u(2*i,2*i);
   }   

   TensorHao<complex<double>, 2> U_Hao(L2,L2),V_Hao(L2,L2);U_Hao=complex <double> (0.0,0.0);V_Hao=complex <double> (0.0,0.0);
   TensorHao<complex<double>, 2> tempMatrix(L2,L2);tempMatrix=complex <double> (0.0,0.0);
   //V_Hao=conj(T_Hao)*lambda_v*transconj(T_Hao);
   BL_NAME(gmm)(conj(T_Hao), lambda_v, tempMatrix ); 
   BL_NAME(gmm)(tempMatrix, T_Hao, V_Hao ,'N','C'); 

   tempMatrix=complex <double> (0.0,0.0);
   //U_Hao=T_Hao*lambda_u*transconj(T_Hao);
   BL_NAME(gmm)(T_Hao, lambda_u, tempMatrix ); 
   BL_NAME(gmm)(tempMatrix, T_Hao, U_Hao ,'N','C'); 

   //get HFB greenMatrix
   TensorHao< complex<double>, 2 > greenMatrix(L2,L2);greenMatrix=complex<double>(0.0,0.0);
   TensorHao< complex<double>, 2 > cmcmMatrix(L2,L2);cmcmMatrix=complex<double>(0.0,0.0);
   TensorHao< complex<double>, 2 > cpcpMatrix(L2,L2);cpcpMatrix=complex<double>(0.0,0.0);
   complex<double> Kenergy(0,0), Venergy(0,0);

   //greenMatrix=(V_Hao^*V_Hao^T)^T
   BL_NAME(gmm)(conj(V_Hao), V_Hao, greenMatrix, 'N','T'); greenMatrix=trans(greenMatrix);

   //cmcmMatrix=(V_Hao^*U_Hao^T)^T
   BL_NAME(gmm)(conj(V_Hao), U_Hao, cmcmMatrix, 'N','T'); cmcmMatrix=trans(cmcmMatrix);

   //cmcmMatrix=-U_Hao^*V_Hao^T
   BL_NAME(gmm)(conj(U_Hao), V_Hao, cpcpMatrix, 'N','T'); cpcpMatrix=complex<double>(-1.0,0.0)*cpcpMatrix;

   //Add K
   const TensorHao< complex<double>, 2 > &K = model->getK();
   for(int i = 0; i < L2; ++i)
   {
       for(int j = 0; j < L2; ++j)
       {
          Kenergy += K(j,i) * greenMatrix(j,i);
       }
   }

   //Add U
   const TensorHao< double, 1> &U = model->getU();
   for(int i = 0; i < L; ++i)
   {
       Venergy += U(i) * ( greenMatrix(i,i)*greenMatrix(i+L,i+L) - greenMatrix(i, i+L)*greenMatrix(i+L, i) );
       Venergy += U(i) * (-1.0* cpcpMatrix(i,i+L)*cmcmMatrix(i,i+L) );
   } 
   double energy = real( Kenergy + Venergy );
   return energy;
}

void GpBCS_VM::get_eiegns(MatrixXcd& cicj_global,MatrixXcd& eigen_values,MatrixXcd& eigen_vectors)
{
    ComplexEigenSolver<MatrixXcd> es( cicj_global ); 
    VectorXcd eigen_values_vector=es.eigenvalues();
       
    eigen_values=MatrixXcd::Zero(2*Nsite,2*Nsite);
    for(int i=1-1;i<=2*Nsite-1;i++)eigen_values(i,i) =  eigen_values_vector(i);
    eigen_vectors = es.eigenvectors();
}  

//stuff to Pure BCS sampling
double GpBCS_VM::get_pure_abEnergy_icf()
{
      MatrixXcd M_jk=MatrixXcd::Zero(Nsite,Nsite);

      for(int j=1-1;j<=Nsite-1;j++){
         for(int k=1-1;k<=Nsite-1;k++){ 
            complex<double> sum_temp=0.0;
            for(int i=1-1;i<=Nsite-1;i++)sum_temp += (T_up(i,j)*T_dn(i,j)*conj(T_up(i,k))*conj(T_dn(i,k)));
            M_jk(j,k) = sum_temp;
         }
      } 

      double energy=0;
      abEnergyRealLimit=0;

      for(int j=1-1;j<=Nsite-1;j++){
         for(int k=1-1;k<=Nsite-1;k++){ 
            energy += real( (C_kl(j,k)/Tot_mark_value)*conj( dia(2*k,2*k+1)/dia(2*j,2*j+1) )*M_jk(j,k) );
            abEnergyRealLimit += -1.0*abs( (C_kl(j,k)/Tot_mark_value)*conj( dia(2*k,2*k+1)/dia(2*j,2*j+1) )*M_jk(j,k) );
         }
      }

      return energy;
}

void GpBCS_VM::get_sampled_C_kl()
{
   int Tot_num=Tot_sample_num;
   Tot_mark_value=0.0;
   C_kl=MatrixXcd::Zero(Nsite,Nsite);

   VectorXi mark=sample_initial_mark();
   double mark_value=sqrt(get_mark_value(mark));
   Tot_mark_value += mark_value;
   update_C_kl(C_kl, mark, mark_value);

   for(int i=2-1;i<=Tot_num-1;i++){
      //std::cout <<" mark check: "<< mark << std::endl;
      sample_update_mark(mark);
      double mark_value=sqrt(get_mark_value(mark));
      Tot_mark_value += mark_value;
      update_C_kl(C_kl, mark, mark_value);
   }
} 

void GpBCS_VM::update_C_kl(MatrixXcd& C_kl, VectorXi& mark, double mark_value)
{
   for(int k=1-1;k<=Nsite-1;k++){
      for(int l=1-1;l<=Nsite-1;l++){
         if(mark(k) == 1 && mark(l) == 0){
           C_kl(k,l) += mark_value;
         }
      }
   }
}







