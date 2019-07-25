#include "../include/vmpBCS.h"

using namespace tensor_hao;
using namespace Eigen;
using namespace LBFGSpp;
using namespace std;

void GpBCS_VM::setMethod(AfqmcConstraintPathMethod* method_)
{
    method=method_;
} 

void GpBCS_VM::setModel(Model* model_, string initialSCPhiTFlag_)
{
    model=model_;
    initialSCPhiTFlag=initialSCPhiTFlag_;
    Nsite=model->getL();
    Ntot=model->getN();
    Nspin=Ntot/2;
    Tot_sample_num=5000;

    if( initialSCPhiTFlag == "setFromDensity_VMGpBCS"){
       N_vm=Nsite;
    }else if( initialSCPhiTFlag == "setFromDensity_VMGpBCS_withGHF_orbital"){
       N_vm=2*Nsite;
    }else if( initialSCPhiTFlag == "setFromDensity_VMGpBCS_withInput_orbital"){
       N_vm=2*Nsite;
    }else if( initialSCPhiTFlag == "setFromGHF"){
       N_vm=2*Nsite;
    }else if( initialSCPhiTFlag == "setFromGHF_readOccupancy"){
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
 
    for(int i=1-1;i<=2*Nsite-1;i++){
        target_dia(i,i)=eigenvalues(i,i); 
    }

    //dia
    for(int i=1-1; i<=Nsite-1; i++){
        complex <double> temp=(target_dia(2*i,2*i)+target_dia(2*i+1,2*i+1))/2.0;
        dia(2*i,2*i+1)=sqrt(temp/(complex <double>(1,0)-temp));
        dia(2*i+1,2*i)=complex<double>(-1.0,0.0)*dia(2*i,2*i+1);
    }
    savedDia=dia;

    updateDia(v);
    
    //T
    if( initialSCPhiTFlag == "setFromDensity_VMGpBCS"){
       for(int i=1-1;i<=2*Nsite-1;i++){
       for(int j=1-1;j<=2*Nsite-1;j++){
          T(i,j)=eigenvectors(i,j);
       }
       }
    }else if( initialSCPhiTFlag == "setFromDensity_VMGpBCS_withGHF_orbital"){
       
       Ghf ghf;
       TensorHao<complex<double>, 2> T_hao=ghf.run();      //icf: we have the error of MPI, need to discuss with Hao

       for(int i=1-1;i<=2*Nsite-1;i++){
       for(int j=1-1;j<=2*Nsite-1;j++){
          T(i,j)=T_hao(i,j);
       }
       }
    }else if( initialSCPhiTFlag == "setFromDensity_VMGpBCS_withInput_orbital"){ 
       
       TensorHao<complex<double>, 2> T_hao;   
       T_hao.read("pBCS_PhiT.T.dat");  

       for(int i=1-1;i<=2*Nsite-1;i++){
       for(int j=1-1;j<=2*Nsite-1;j++){
          T(i,j)=T_hao(i,j);
       }
       }
    }else if( initialSCPhiTFlag == "setFromGHF"){
       
       Ghf ghf;
       TensorHao<complex<double>, 2> T_hao=ghf.run();      

       for(int i=1-1;i<=2*Nsite-1;i++){
       for(int j=1-1;j<=2*Nsite-1;j++){
          T(i,j)=T_hao(i,j);
       }
       }
    }else if( initialSCPhiTFlag == "setFromGHF_readOccupancy"){
       
       Ghf ghf;
       TensorHao<complex<double>, 2> T_hao=ghf.run();      

       for(int i=1-1;i<=2*Nsite-1;i++){
       for(int j=1-1;j<=2*Nsite-1;j++){
          T(i,j)=T_hao(i,j);
       }
       }
    }else{
      cout<<"initialSCPhiTFlag Error2: "<<initialSCPhiTFlag<<endl;
      exit(1);
    }

    //Rotation
    TensorHao<complex<double>, 2> Rotation(2*Nsite,2*Nsite); Rotation = complex <double> (0.0,0.0);  //Attention: only work for "sawtooth" coupling
    TensorHao<complex<double>, 2> T_hao(2*Nsite,2*Nsite); T_hao = complex <double> (0.0,0.0);
    TensorHao<complex<double>, 2> T_R(2*Nsite,2*Nsite); T_R = complex <double> (0.0,0.0);       

    for(int i=1-1;i<=2*Nsite-1;i++){
    for(int j=1-1;j<=2*Nsite-1;j++){
       T_hao(i,j)=T(i,j);
    }
    }
    for(int i=1-1;i<=Nsite-1;i++){
       Rotation(2*i,2*i)=method->rotation00;
       Rotation(2*i,2*i+1)=method->rotation01;
       Rotation(2*i+1,2*i)=method->rotation10;
       Rotation(2*i+1,2*i+1)=method->rotation11;
    }
    BL_NAME(gmm)(Rotation, T_hao , T_R); 

    for(int i=1-1;i<=2*Nsite-1;i++){
    for(int j=1-1;j<=2*Nsite-1;j++){
       T(i,j)=T_R(i,j);
    }
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
      getEnergyFast();
   }else if( initialSCPhiTFlag == "setFromDensity_VMGpBCS_withGHF_orbital"){
      //get density matrix error from GpBCS (<GpBCS|CiCj|GpBCS>/<GpBCS|GpBCS>)
      getDensityMatrixDiatance();
   }else if( initialSCPhiTFlag == "setFromDensity_VMGpBCS_withInput_orbital"){
      //get density matrix error from GpBCS (<GpBCS|CiCj|GpBCS>/<GpBCS|GpBCS>)
      getDensityMatrixDiatance();
   }else if( initialSCPhiTFlag == "setFromGHF"){
      //get energy from GpBCS (<GpBCS|CiCj|GpBCS>/<GpBCS|GpBCS>)
      getEnergyFast();
   }else if( initialSCPhiTFlag == "setFromGHF_readOccupancy"){
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
   }else if( initialSCPhiTFlag == "setFromDensity_VMGpBCS_withGHF_orbital"){
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
   }else if( initialSCPhiTFlag == "setFromDensity_VMGpBCS_withInput_orbital"){  
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
   }else if( initialSCPhiTFlag == "setFromGHF_readOccupancy"){
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

void GpBCS_VM::get_eiegns(MatrixXcd& cicj_global,MatrixXcd& eigen_values,MatrixXcd& eigen_vectors)
{
    ComplexEigenSolver<MatrixXcd> es( cicj_global ); 
    VectorXcd eigen_values_vector=es.eigenvalues();
       
    eigen_values=MatrixXcd::Zero(2*Nsite,2*Nsite);
    for(int i=1-1;i<=2*Nsite-1;i++)eigen_values(i,i) =  eigen_values_vector(i);
    eigen_vectors = es.eigenvectors();
}  
