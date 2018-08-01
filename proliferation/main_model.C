#include <iostream>
#include <stdio.h>
#include <fstream>
#include <string.h>
#include <sstream>
#include <math.h>
#include <vector>

#include "main_model.h"

using namespace std;

void RungeKutta(std::vector<double>& QoI,const double time_step,vector<double> Parameters);
void fmodel(std::vector<double> qoi,std::vector<double>& auxiliar,vector<double> Parameters);

void main_code(vector<double>& TUMOR_MASS,vector<double> Parameters,double tumor_mean,vector<double> SaveTime){
  
  /////----- Read argumments -----/////
  const unsigned int n_times     = SaveTime.size();
  const double time_step         = 0.01;
  const double l_day             = SaveTime[n_times-1];
  const double init_time         = SaveTime[0];
  const unsigned int n_timesteps = (l_day*100.0-init_time*100.0)/(time_step*100.0);
  std::vector<unsigned int> save_iter(n_times,0);
  for(unsigned int i = 0; i < n_times; i++)
    save_iter[i] = (SaveTime[i]*100.0-init_time*100.0)/(time_step*100.0);

  const unsigned int number_equations = 1;
  std::vector<double> QoI(number_equations,0.0);
  
  QoI[0] = tumor_mean;
  TUMOR_MASS.push_back(QoI[0]);
  /////----- Defining parameters -----/////
  unsigned int t_step = 0;

  /////----- Solving the problem -----/////
  do
    {
      /////----- Prepare time step-----/////
      t_step++;
      RungeKutta(QoI,time_step,Parameters);
      
      /////----- Post-processing -----/////
      int need_to_print = 0;
      for(unsigned int i = 0; i < n_times; i++)
	if(save_iter[i] == t_step)
	  need_to_print = 1;
      if(need_to_print)
	TUMOR_MASS.push_back(QoI[0]);
    }
  while(t_step<n_timesteps);
}

void RungeKutta(std::vector<double>& QoI,const double time_step,vector<double> Parameters){

  const unsigned int N = QoI.size();
  std::vector<double> k(4*N,0.0);
  std::vector<double> qoi(N,0.0);
  std::vector<double> auxiliar(N,0.0);
  
  for(unsigned int i = 0; i < N; i++){qoi[i] = QoI[i];}
  fmodel(qoi,auxiliar,Parameters);
  for(unsigned int i = 0; i < N; i++){k[i*4+0] = time_step*auxiliar[i];}
  
  for(unsigned int i = 0; i < N; i++){qoi[i] = QoI[i]+(k[i*4+0]/2.);}
  fmodel(qoi,auxiliar,Parameters);
  for(unsigned int i = 0; i < N; i++){k[i*4+1] = time_step*auxiliar[i];}
 
  for(unsigned int i = 0; i < N; i++){qoi[i] = QoI[i]+(k[i*4+1]/2.);}
  fmodel(qoi,auxiliar,Parameters);
  for(unsigned int i = 0; i < N; i++){k[i*4+2] = time_step*auxiliar[i];}

  for(unsigned int i = 0; i < N; i++){qoi[i] = QoI[i]+k[i*4+2];}
  fmodel(qoi,auxiliar,Parameters);
  for(unsigned int i = 0; i < N; i++){k[i*4+3] = time_step*auxiliar[i];}
  
  for(unsigned int i = 0; i < N; i++){QoI[i] = QoI[i]+(1./6.)*(k[i*4+0]+2*(k[i*4+1]+k[i*4+2])+k[i*4+3]);}
  
}

void fmodel(std::vector<double> qoi,std::vector<double>& auxiliar,vector<double> Parameters){
  double prol = Parameters[0];
  double supp = Parameters[1];
  double apop = Parameters[2];
  auxiliar[0]=prol*qoi[0]*(1.0-qoi[0]/supp)-apop*qoi[0];
}
