//
// Created by boruoshihao on 4/16/17.
// Modefied by icf on 05/01/2019
//
#include "../include/afqmcConstraintPath.h"
#include "../include/vmpBCS.h"

using namespace std;
using namespace tensor_hao;

using namespace Eigen;
using namespace LBFGSpp;

////////////////////////////////////////////////////////////////////////////////
//need static function in VM
////////////////////////////////////////////////////////////////////////////////

GpBCS_VM gpBCS_VM;  //icf: because of muti-define, we have to do so. Discuss with Hao to solve it. !!!!!!!!!!!!!!!!!!!!!!! 

double VM(VectorXd& v,VectorXd& vgrad)
{
    double step=0.001;
    
    v(0)=0.0; //icf: phase lock
    for(int i=2-1;i<=gpBCS_VM.N_vm-1;i++){
       if(real(v(i)) > 10.0)v(i)=10.0;
       if(real(v(i)) < -10.0)v(i)=-10.0;
    }

    VectorXcd v_d=VectorXcd::Zero(gpBCS_VM.N_vm);
    for(int i=1-1;i<=gpBCS_VM.N_vm-1;i++)v_d(i)=v(i);


    double energy_in_VM=gpBCS_VM.e_v_update(v_d);
    cout<<"energy_in_VM: "<<energy_in_VM<<endl;

    for(int i=1-1;i<=gpBCS_VM.N_vm-1;i++)v(i)=real(v_d(i));


    VectorXcd v_change1=v_d;
    VectorXcd v_change2=v_d;
    
    vgrad( 0  )=0.0;
    for(int k=2-1;k<=gpBCS_VM.N_vm-1;k++){    //icf: phase lock
       v_change1=v_d;
       v_change2=v_d;
       v_change1( k )=v_change1( k  )+step;
       v_change2( k )=v_change2( k  )-step;
       vgrad( k  )=(gpBCS_VM.e_v_update(v_change1)-gpBCS_VM.e_v_update(v_change2))/(2.0*step);
    }
    return energy_in_VM;
}

void vmpBCS(const TensorHao<complex<double>, 2> &tempGreen,TensorHao<complex<double>, 2> &tempOrbital,TensorHao<complex<double>, 2> &tempAnalytic)
{
    cout<<"get into VM"<<endl;
    gpBCS_VM.outfile_log.open("LOG.dat");
    gpBCS_VM.startTime=clock();
    //start log

    MatrixXcd dm=MatrixXcd::Zero(2*gpBCS_VM.Nsite,2*gpBCS_VM.Nsite);
    MatrixXcd dm_dia=MatrixXcd::Zero(2*gpBCS_VM.Nsite,2*gpBCS_VM.Nsite);
    MatrixXcd dm_T=MatrixXcd::Zero(2*gpBCS_VM.Nsite,2*gpBCS_VM.Nsite);
    for(int i=1-1;i<=2*gpBCS_VM.Nsite-1;i++){
       for(int j=1-1;j<=2*gpBCS_VM.Nsite-1;j++){
          dm(i,j)=tempGreen(i,j);
          dm_dia(i,j)=tempAnalytic(i,j);
          dm_T(i,j)=tempOrbital(i,j);
       }
    }
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //////GpBCS///////////////////////////////////////
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VectorXcd v_d=VectorXcd::Zero(gpBCS_VM.N_vm);
    VectorXd v=VectorXd::Zero(gpBCS_VM.N_vm);
    VectorXd vgrad=VectorXd::Zero(gpBCS_VM.N_vm);

    LBFGSParam<double> param;
    param.max_iterations = 100;     //icf: param used for mixedVM only!!!
    param.epsilon=0.001; 
    //param.delta=0.001;
    //param.past=1;
    
    LBFGSSolver<double> solver(param);

    if( gpBCS_VM.initialSCPhiTFlag == "setFromDensity_VMGpBCS"){
       //for(int i=1-1; i<=gpBCS_VM.N_vm-1; i++)v(i)=0.0;
       for(int i=1-1; i<=gpBCS_VM.N_vm-1; i++)v(i)=2.0*rand()/double(RAND_MAX)-1.0;
    }else if( gpBCS_VM.initialSCPhiTFlag == "setFromGHF"){
       for(int i=1-1; i<=gpBCS_VM.N_vm/2-1; i++){
           v(2*i)=0.0;
           if(i<=gpBCS_VM.Nsite-gpBCS_VM.Nspin-1){
              v(2*i+1)=0.0;  
           }else{
              v(2*i+1)=1.0; 
           }        
       }
    }else{
       cout<<"initialSCPhiTFlag Error: "<<gpBCS_VM.initialSCPhiTFlag<<endl;
       exit(1);
    }

    for(int i=1-1;i<=gpBCS_VM.N_vm-1;i++)v_d(i)=v(i);

    gpBCS_VM.set(dm, dm_dia, dm_T, v_d);

    if( gpBCS_VM.initialSCPhiTFlag == "setFromGHF" ){
    }else if( gpBCS_VM.initialSCPhiTFlag == "setFromDensity_VMGpBCS"){

       gpBCS_VM.get_sampled_C_kl(); //ATTENTION: need to set C_jk first (before VM) for "getEnergyFastICF" and comment out for efficiency if you use other energyEstimator
    
       double fw=0;
       int niter = solver.minimize(VM, v, fw);
       std::cout << niter << " iterations "<< fw << std::endl;

       cout<<"Be carefull: we are using abEnergy_icf, which means pBCS T_up and T_dn are used --> make sure you got correct T_up and T_dn for GpBCS"<<endl;
       cout<<"the limit we may get: "<<gpBCS_VM.abEnergyRealLimit<<endl;
    }
    std::cout << "v: \n"<< v << std::endl;

    for(int i=1-1;i<=gpBCS_VM.N_vm-1;i++)v_d(i)=v(i);
    
    gpBCS_VM.e_v_update(v_d);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //////End GpBCS///////////////////////////////////////
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    tempOrbital.resize(2*gpBCS_VM.Nsite,2*gpBCS_VM.Nsite);tempAnalytic.resize(2*gpBCS_VM.Nsite,2*gpBCS_VM.Nsite);
    for(int i=1-1;i<=2*gpBCS_VM.Nsite-1;i++){
       for(int j=1-1;j<=2*gpBCS_VM.Nsite-1;j++){
          tempOrbital(i,j)=gpBCS_VM.T(i,j);
          tempAnalytic(i,j)=gpBCS_VM.dia(i,j);
       }
    }

    //end log
    gpBCS_VM.outfile_log.close();    
}
////////////////////////////////////////////////////////////////////////////////
//need static function in VM
////////////////////////////////////////////////////////////////////////////////


void AfqmcConstraintPath::initialPhiT()
{
    if(method.initialPhiTFlag == "setFromModel") //icf: we could do betterrrr here
    {
        if( MPIRank()==0 ){
           size_t L=2*model.getL(); 
           size_t N=model.getN();

           TensorHao<complex<double>, 2> tempOrbital_flip = model.getKEigenVector();
           //cout<<model.getKEigenValue()<<endl;exit(1);
           TensorHao<complex<double>, 2> tempOrbital(L,L);
           TensorHao<complex<double>, 2> tempAnalytic(L,L);
           TensorHao<complex<double>, 2> tempMatrix(L,L);
           TensorHao<complex<double>, 2> tempMatrix2(L,L);

           for(size_t i=1-1; i<=L/2-N/2-1; i++){
              tempAnalytic(2*i,2*i+1)=complex<double>(0.0,0.0);
              tempAnalytic(2*i+1,2*i)=complex<double>(-1.0,0.0)*tempAnalytic(2*i,2*i+1);
           }
           for(size_t i=L/2-N/2; i<=L/2-1; i++){
              tempAnalytic(2*i,2*i+1)=complex<double>(1.0,0.0);
              tempAnalytic(2*i+1,2*i)=complex<double>(-1.0,0.0)*tempAnalytic(2*i,2*i+1);
           }

           for(size_t j=1-1; j<=L-1; j++){ 
           for(size_t i=1-1; i<=L-1; i++){ 
              tempOrbital(i,j)=tempOrbital_flip(i,L-1-j);
           }
           }

           BL_NAME(gmm)(tempOrbital, tempAnalytic, tempMatrix ); 
           BL_NAME(gmm)(tempMatrix, tempOrbital, tempMatrix2 ,'N','T'); 

           phiT.wfRef()=tempMatrix2;
           phiT.logwRef()=complex<double>(0.0,0.0);
           phiT.orbital=tempOrbital;
           phiT.occupancy=tempAnalytic;
        } 

        if( MPIRank()==0 ){
           phiT.write("pBCS_PhiT.Wf.dat0");
           phiT.writeOrbital("pBCS_PhiT.T.dat0");
           phiT.writeOccupancy("pBCS_PhiT.Gdia.dat0");
        }
    }
    else if(method.initialPhiTFlag == "setRandomly")
    {
        //if( MPIRank()==0 ) fillWalkerRandomly(phiT, model);
        //if( MPIRank()==0 ) phiT.write("phiT.dat");
        //MPIBcast(phiT);
        cout<<"Error!!! wo don't do setRandomly!"<<endl;
        exit(1);
    }
    else if(method.initialPhiTFlag == "readFromFile")
    {
        if( MPIRank()==0 ){
           size_t L=2*model.getL(); 
           size_t N=model.getN();

           TensorHao<complex<double>, 2> tempMatrix(L,L);
           TensorHao<complex<double>, 2> tempMatrix2(L,L);

           phiT.read("pBCS_PhiT.Wf.dat");
           phiT.readOrbital("pBCS_PhiT.T.dat");
           phiT.readOccupancy("pBCS_PhiT.Gdia.dat");

           BL_NAME(gmm)(phiT.orbital, phiT.occupancy, tempMatrix ); 
           BL_NAME(gmm)(tempMatrix, phiT.orbital, tempMatrix2 ,'N','T');

           phiT.wfRef()=tempMatrix2;
           phiT.logwRef()=complex<double>(0.0,0.0);
        }
    }
    else
    {
        cout<<"Error!!! Do not recognize initialPhiTFlag!"<<endl;
        exit(1);
    }
    MPIBcast(phiT);
}

void AfqmcConstraintPath::initialWalker()
{
    walker.resize(method.walkerSizePerThread);
    walkerIsAlive.resize(method.walkerSizePerThread);

    size_t N=model.getN();

    if(method.initialWalkerFlag == "setFromModel")
    {
        if( MPIRank()==0 ) fillWalkerFromModel(walker[0], model);

        if( MPIRank()==0){
           WalkerWalkerOperation walkerWalkerOperation;
           walkerWalkerOperation.set(phiT, walker[0]);  
           complex <double> overlap = exp( walkerWalkerOperation.returnLogOverlap() );   
           complex <double> phase = ( walkerWalkerOperation.returnLogOverlap() ).imag();

           cout<<"We get ovelap: "<<overlap<<endl;
           cout<<"We get phase: "<<phase<<endl;
 
           phiT.wfRef()=exp(complex <double> (0,2.0/double(N))*phase)*phiT.wfRef();     //icf: pf(phi^T M^* phi=N*N matrix)
           phiT.occupancy=exp(complex <double> (0,2.0/double(N))*phase)*phiT.occupancy;

           walkerWalkerOperation.set(phiT, walker[0]);  
           overlap = exp( walkerWalkerOperation.returnLogOverlap() );  
           cout<<"We modify the phiT: "<<overlap<<endl;

        }
        MPIBcast(phiT);

        if( MPIRank()==0 ) walker[0].write("phi.dat");
        MPIBcast(walker[0]);

        for(int i = 0; i < method.walkerSizePerThread; ++i)
        {
            walker[i] = walker[0];
            walkerIsAlive[i] = true;
        }

    }
    else if(method.initialWalkerFlag == "setRandomly")
    {
        if( MPIRank()==0 ) fillWalkerRandomly(walker[0], model);
        if( MPIRank()==0 ) walker[0].write("phi.dat");
        MPIBcast(walker[0]);

        for(int i = 0; i < method.walkerSizePerThread; ++i)
        {
            walker[i] = walker[0];
            walkerIsAlive[i] = true;
        }
    }
    else if(method.initialWalkerFlag == "sampleFromPhiT")
    {
        if( MPIRank()==0 ){
           setWalkerFromPhiT(walker[0], phiT, 2*model.getL(), model.getN());  

           WalkerWalkerOperation walkerWalkerOperation;
           walkerWalkerOperation.set(phiT, walker[0]);  
           complex <double> overlap = exp( walkerWalkerOperation.returnLogOverlap() );   
           complex <double> phase = ( walkerWalkerOperation.returnLogOverlap() ).imag();

           cout<<"We get ovelap: "<<overlap<<endl;
           cout<<"We get phase: "<<phase<<endl;
 
           phiT.wfRef()=exp(complex <double> (0,2.0/double(N))*phase)*phiT.wfRef();     //icf: pf(phi^T M^* phi=N*N matrix)
           phiT.occupancy=exp(complex <double> (0,2.0/double(N))*phase)*phiT.occupancy;

           walkerWalkerOperation.set(phiT, walker[0]);  
           overlap = exp( walkerWalkerOperation.returnLogOverlap() );  
           cout<<"We modify the phiT: "<<overlap<<endl;
        }
        MPIBcast(phiT);

        if( MPIRank()==0 ) walker[0].write("phi.dat");
        MPIBcast(walker[0]);

        for(int i = 0; i < method.walkerSizePerThread; ++i)
        {
            walker[i] = walker[0];
            walkerIsAlive[i] = true;
        }
    }
    else if(method.initialWalkerFlag == "readFromFile")
    {
        if( MPIRank()==0 ) walker[0].read("phi.dat");
        MPIBcast(walker[0]);

        for(int i = 0; i < method.walkerSizePerThread; ++i)
        {
            walker[i] = walker[0];
            walkerIsAlive[i] = true;
        }
    }
    else if(method.initialWalkerFlag == "readAllWalkers")
    {
        string filename;
        int baseNumber = MPIRank() * method.walkerSizePerThread;
        for(int i = 0; i < method.walkerSizePerThread; ++i)
        {
            filename = "./walkers/phi_" + to_string(i+baseNumber) +".dat";
            if( checkFile(filename) ) { walker[i].read(filename); walkerIsAlive[i]=true; }
            else { fillWalkerRandomly(walker[i], model);  walkerIsAlive[i]=false; }
        }
    }
    else
    {
        cout<<"Error!!! Do not recognize initialWalkerFlag!"<<endl;
        exit(1);
    }

    initialMgsAndPopControl();

}

void AfqmcConstraintPath::initialSCPhiT()
{
  if( MPIRank()==0 ){
     size_t L=2*model.getL();size_t N=model.getN();
     TensorHao<complex<double>, 2> tempGreen(L,L); tempGreen=complex <double> (0.0,0.0);
     TensorHao<complex<double>, 2> tempGreenUp(L/2,L/2); tempGreen=complex <double> (0.0,0.0);
     TensorHao<complex<double>, 2> tempGreenDn(L/2,L/2); tempGreen=complex <double> (0.0,0.0);
     TensorHao<double, 1> tempDia(L),tempDiaUp(L/2),tempDiaDn(L/2); tempDiaUp=0.0;tempDiaDn=0.0;
     TensorHao<complex<double>, 2> tempOrbital(L,L);tempOrbital=complex <double> (0.0,0.0);
     TensorHao<complex<double>, 2> tempAnalytic(L,L);tempAnalytic(L,L)=complex <double> (0.0,0.0);

     for(size_t i=1-1; i<=L/2-1; i++){
     for(size_t j=1-1; j<=L/2-1; j++){
        tempGreen(i,j)=greenMatrixAfterCPMC(i,j);
        tempGreenUp(i,j)=tempGreen(i,j);
        tempGreen(i+L/2,j+L/2)=greenMatrixAfterCPMC(i+L/2,j+L/2);
        tempGreenDn(i,j)=tempGreen(i+L/2,j+L/2);
     }
     }
     tempGreen=(tempGreen+conjtrans(tempGreen))/(complex <double> (2,0));
     tempGreenUp=(tempGreenUp+conjtrans(tempGreenUp))/(complex <double> (2,0));
     tempGreenDn=(tempGreenDn+conjtrans(tempGreenDn))/(complex <double> (2,0));

     checkHermitian(tempGreenUp);
     checkHermitian(tempGreenDn);
     //U=U*D*V
     BL_NAME(eigen)(tempGreenUp,tempDiaUp);
     BL_NAME(eigen)(tempGreenDn,tempDiaDn);

     for(size_t i=1-1; i<=L/2-1; i++){
     for(size_t j=1-1; j<=L/2-1; j++){
        tempOrbital(j,2*i)=tempGreenUp(j,i);
        tempDia(2*i)=tempDiaUp(i);
        tempOrbital(j+L/2,2*i+1)=tempGreenDn(j,i);
        tempDia(2*i+1)=tempDiaDn(i);
     }
     }

     cout<<"The eigen of GreenMatrix: "<<endl;
     cout<<tempDia<<endl; 

    if(method.initialSCPhiTFlag == "setFromDensity_Analytical") //icf: we could do betterrrr here
    {
   
           for(size_t i=1-1; i<=L/2-1; i++){
              complex <double> temp=complex <double>((tempDia(2*i)+tempDia(2*i+1))/2.0);
              tempAnalytic(2*i,2*i+1)=sqrt(temp/(complex <double>(1,0)-temp));
              tempAnalytic(2*i+1,2*i)=complex<double>(-1.0,0.0)*tempAnalytic(2*i,2*i+1);
           }

           TensorHao<complex<double>, 2> tempMatrix(L,L);
           TensorHao<complex<double>, 2> tempMatrix2(L,L);

           BL_NAME(gmm)(tempOrbital, tempAnalytic, tempMatrix ); 
           BL_NAME(gmm)(tempMatrix, tempOrbital, tempMatrix2 ,'N','T'); 

           phiT.wfRef()=tempMatrix2;
           phiT.logwRef()=complex<double>(0.0,0.0);
           phiT.orbital=tempOrbital;
           phiT.occupancy=tempAnalytic;


          ostringstream oss1;
          oss1 << "pBCS_PhiT.Wf.dat0"<<scLooPSteps;

          ostringstream oss2;
          oss2 << "pBCS_PhiT.T.dat0"<<scLooPSteps;

          ostringstream oss3;
          oss3 << "pBCS_PhiT.Gdia.dat0"<<scLooPSteps;

           phiT.write(oss1.str());
           phiT.writeOrbital(oss2.str());
           phiT.writeOccupancy(oss3.str());
    }
    else if(method.initialSCPhiTFlag == "setFromDensity_DET") 
    {

           for(size_t i=L/2-N/2; i<=L/2-1; i++){
              cout<<"set "<<tempDia(2*i)<<" and "<<tempDia(2*i+1)<<" to be "<<1.0<<endl;
              tempAnalytic(2*i,2*i+1)=1.0;
              tempAnalytic(2*i+1,2*i)=complex<double>(-1.0,0.0)*tempAnalytic(2*i,2*i+1);
           }

           TensorHao<complex<double>, 2> tempMatrix(L,L);
           TensorHao<complex<double>, 2> tempMatrix2(L,L);

           BL_NAME(gmm)(tempOrbital, tempAnalytic, tempMatrix ); 
           BL_NAME(gmm)(tempMatrix, tempOrbital, tempMatrix2 ,'N','T'); 

           phiT.wfRef()=tempMatrix2;
           phiT.logwRef()=complex<double>(0.0,0.0);
           phiT.orbital=tempOrbital;
           phiT.occupancy=tempAnalytic;

          ostringstream oss1;
          oss1 << "pBCS_PhiT.Wf.dat0"<<scLooPSteps;

          ostringstream oss2;
          oss2 << "pBCS_PhiT.T.dat0"<<scLooPSteps;

          ostringstream oss3;
          oss3 << "pBCS_PhiT.Gdia.dat0"<<scLooPSteps;

           phiT.write(oss1.str());
           phiT.writeOrbital(oss2.str());
           phiT.writeOccupancy(oss3.str());
    }
    else if(method.initialSCPhiTFlag == "setFromDensity_VMGpBCS") //icf: we could do betterrrr here
    {

           for(size_t i=1-1; i<=L/2-1; i++){
              complex <double> temp=complex <double>((tempDia(2*i)+tempDia(2*i+1))/2.0);
              tempAnalytic(2*i,2*i+1)=sqrt(temp/(complex <double>(1,0)-temp));
              tempAnalytic(2*i+1,2*i)=complex<double>(-1.0,0.0)*tempAnalytic(2*i,2*i+1);
           }

           gpBCS_VM.setModel(&model,method.initialSCPhiTFlag);
           vmpBCS(tempGreen,tempOrbital,tempAnalytic);

           TensorHao<complex<double>, 2> tempMatrix(L,L);
           TensorHao<complex<double>, 2> tempMatrix2(L,L);

           BL_NAME(gmm)(tempOrbital, tempAnalytic, tempMatrix ); 
           BL_NAME(gmm)(tempMatrix, tempOrbital, tempMatrix2 ,'N','T'); 

           phiT.wfRef()=tempMatrix2;
           phiT.logwRef()=complex<double>(0.0,0.0);
           phiT.orbital=tempOrbital;
           phiT.occupancy=tempAnalytic;


          ostringstream oss1;
          oss1 << "pBCS_PhiT.Wf.dat"<<scLooPSteps;

          ostringstream oss2;
          oss2 << "pBCS_PhiT.T.dat"<<scLooPSteps;

          ostringstream oss3;
          oss3 << "pBCS_PhiT.Gdia.dat"<<scLooPSteps;

           phiT.write(oss1.str());
           phiT.writeOrbital(oss2.str());
           phiT.writeOccupancy(oss3.str());

    }
    else if(method.initialSCPhiTFlag == "setFromGHF") 
    {
           size_t L=2*model.getL();
           
           gpBCS_VM.setModel(&model,method.initialSCPhiTFlag);
           vmpBCS(tempGreen,tempOrbital,tempAnalytic);

           TensorHao<complex<double>, 2> tempMatrix(L,L);
           TensorHao<complex<double>, 2> tempMatrix2(L,L);

           BL_NAME(gmm)(tempOrbital, tempAnalytic, tempMatrix ); 
           BL_NAME(gmm)(tempMatrix, tempOrbital, tempMatrix2 ,'N','T'); 

           phiT.wfRef()=tempMatrix2;
           phiT.logwRef()=complex<double>(0.0,0.0);
           phiT.orbital=tempOrbital;
           phiT.occupancy=tempAnalytic;


          ostringstream oss1;
          oss1 << "pBCS_PhiT.Wf.dat"<<scLooPSteps;

          ostringstream oss2;
          oss2 << "pBCS_PhiT.T.dat"<<scLooPSteps;

          ostringstream oss3;
          oss3 << "pBCS_PhiT.Gdia.dat"<<scLooPSteps;

           phiT.write(oss1.str());
           phiT.writeOrbital(oss2.str());
           phiT.writeOccupancy(oss3.str());

    }
    else if(method.initialSCPhiTFlag == "setRandomly")
    {
        //if( MPIRank()==0 ) fillWalkerRandomly(phiT, model);
        //if( MPIRank()==0 ) phiT.write("phiT.dat");
        //MPIBcast(phiT);
        cout<<"Error!!! wo don't do setRandomly!"<<endl;
        exit(1);
    }
    else
    {
        cout<<"Error!!! Do not recognize initialPhiTSCFlag!"<<endl;
        exit(1);
    }
  } 
  MPIBcast(phiT);
}

void AfqmcConstraintPath::writeWalkers()
{
    MPIBarrier();
    if( MPIRank() == 0 )
    {
        int flag;
        flag = system("rm -rf walkers");   if(flag != 0) cout<<"WARNING!!! system command does not exit properly!"<<endl;
        flag = system("mkdir -p walkers"); if(flag != 0) cout<<"WARNING!!! system command does not exit properly!"<<endl;
    }
    MPIBarrier();

    string filename;
    int baseNumber = MPIRank() * method.walkerSizePerThread;
    for(int i = 0; i < method.walkerSizePerThread; ++i)
    {
        if( walkerIsAlive[i] )
        {
            filename = "./walkers/phi_" + to_string(i+baseNumber) +".dat";
            walker[i].write(filename);
        }
    }
}

void AfqmcConstraintPath::initialMgsAndPopControl()
{

    modifyGM();
    popControl();
}



