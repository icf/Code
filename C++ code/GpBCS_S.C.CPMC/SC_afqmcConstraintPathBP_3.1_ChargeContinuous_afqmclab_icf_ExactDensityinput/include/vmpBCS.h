#ifndef VMPBCS
#define VMPBCS

#define _USE_MATH_DEFINES

#include<iostream>
#include<fstream>

//#include <complex.h>
#include <iomanip>

#include <math.h>
#include<cmath>
#include <stdlib.h> 

#include <vector> 
//#include <Eigen/Core>
//#include <Eigen/Dense>
//#include <Eigen/Eigenvalues>
//#include <unsupported/Eigen/MatrixFunctions>

#include <algorithm>

#include <unistd.h>

#include<time.h>

#include "../include/eigen/Eigen/Core"
#include "../include/eigen/Eigen/Dense"
#include "../include/eigen/Eigen/Eigenvalues"
#include "../include/LBFGSpp-master/include/LBFGS.h"
#include "../include/afqmcConstraintPath.h"
#include "afqmclab.h"
#include "../include/ghf.h"

using namespace tensor_hao;
using namespace Eigen;
using namespace LBFGSpp;
using namespace std;

class GpBCS_VM
{
   public: 
   ofstream outfile_log;
   clock_t startTime,endTime;

   int Nsite,Ntot,Nspin,N_vm,Tot_sample_num; 
   double Tot_mark_value;

   MatrixXcd dia,savedDia,T,M,target_dia; 
     
   MatrixXcd target_dm,dm,target_nn,nn;
   double energy;
     
   Model* model;
   AfqmcConstraintPathMethod* method;
   string initialSCPhiTFlag;

   void setModel(Model* model_, string initialSCPhiTFlag_); 
   void setMethod(AfqmcConstraintPathMethod* method_); 
   void set(MatrixXcd& target_dm_, MatrixXcd& eigenvalues, MatrixXcd& eigenvectors, VectorXcd& v);
   //get sampled SD from GpBCS
   VectorXi sample_initial_mark(void);
   double get_mark_value(VectorXi& mark);
   MatrixXcd get_mark_SD(VectorXi& mark);
   void sample_update_mark(VectorXi& mark);
   //get pure energy estimator
   void updateDia(VectorXcd& v); 
   double e_v_update(VectorXcd& v);

   void getEnergyFast();
   void getEnergyFastMixed();
   void getDensityMatrixDiatance();

   double get_mixed_energy(VectorXi& mark);
   MatrixXcd get_mixed_dm(VectorXi& mark, double &Energy);

   //tool
   void get_eiegns(MatrixXcd& cicj_global,MatrixXcd& eigen_values,MatrixXcd& eigen_vectors);

};// end HSystem
#endif //VMPBCS

