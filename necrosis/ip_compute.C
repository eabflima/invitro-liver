#include <cmath>
#include <sys/time.h>

#include <queso/GenericVectorFunction.h>
#include <queso/GaussianVectorRV.h>
#include <queso/GenericVectorRV.h>
#include <queso/ConcatenationSubset.h>
#include <queso/ConcatenatedVectorRV.h>

#include <queso/GenericScalarFunction.h>
#include <queso/GslMatrix.h>
#include <queso/UniformVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/StatisticalForwardProblem.h>

#include <ip_compute.h>
#include <ip_likelihood.h>
#include <fp_qoi.h>

#include "libmesh/libmesh.h"
#include "libmesh/getpot.h"

using namespace libMesh;
using namespace std;

void computeIP(const QUESO::FullEnvironment& env)
{
  //------------------------------------------------------
  //
  //                 CALIBRATION STEP
  //
  //------------------------------------------------------
  cout << endl << endl << "                       Starting Calibration Step" << endl << endl;
  //------------------------------------------------------
  // SIP Step 1 of 6: Instantiate the parameter space
  //------------------------------------------------------
  GetPot input_data("inp_dat.in");

  unsigned int N_par=6;

  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix>
    paramSpaceA(env, "paramA_", 2, NULL);
  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix>
    paramSpaceB(env, "paramB_", 4, NULL);
  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix>
    paramSpace (env, "param_", N_par, NULL);
  
  //------------------------------------------------------
  // SIP Step 2 of 6: Instantiate the parameter domain
  //------------------------------------------------------
  QUESO::GslVector paramMinsA(paramSpaceA.zeroVector());
  paramMinsA[0] = input_data("v_l_n_min",0.10);
  paramMinsA[1] = input_data("v_std_min",10.0);
  QUESO::GslVector paramMaxsA(paramSpaceA.zeroVector());
  paramMaxsA[0] = input_data("v_l_n_max",0.55);
  paramMaxsA[1] = input_data("v_std_max",20000.0);
  QUESO::BoxSubset<QUESO::GslVector,QUESO::GslMatrix>
    paramDomainA("paramA_",paramSpaceA,paramMinsA,paramMaxsA);
  
  QUESO::GslVector paramMinsB(paramSpaceB.zeroVector());
  paramMinsB.cwSet(0);
  QUESO::GslVector paramMaxsB(paramSpaceB.zeroVector());
  paramMaxsB.cwSet(INFINITY);
  QUESO::BoxSubset<QUESO::GslVector,QUESO::GslMatrix>
    paramDomainB("paramB_",paramSpaceB,paramMinsB,paramMaxsB);

  QUESO::ConcatenationSubset<QUESO::GslVector,QUESO::GslMatrix>
    paramDomain("param_",paramSpace,paramDomainA,paramDomainB);

  //------------------------------------------------------
  // SIP Step 3 of 6: Instantiate the likelihood function
  // object to be used by QUESO.
  //------------------------------------------------------

  // Data available in cal_data.dat
  std::ifstream read;
  std::string line;
  double mass,time;
  std::vector<double> v_mass;
  std::vector<double> v_time;
  std::string cal_name;
  cal_name.append(input_data("cal_name", "cal_data.dat"));
  const unsigned int n_samples = input_data("n_samples",4);
  const unsigned int n_cells   = input_data("n_cells",5);
  const unsigned int cell_type = input_data("cell_type",1);
  read.open(cal_name);
  if(!read.is_open())
    std::cout << "Error opening data file.\n";
  unsigned int count = 0;
  while(std::getline(read,line)){
    std::istringstream iss (line);
    count++;
    iss >> time;
    if(count == 1){
       v_time.push_back(time);
    }
    else if(count == n_samples)
      count = 0;
    for(unsigned int p = 1; p <= n_cells; p++){
      iss >> mass;
      if(p == cell_type)
	v_mass.push_back(mass);
    }
  }
  read.close();

  unsigned int Size = v_time.size();

  likelihoodRoutine_DataType likelihoodRoutine_Data;
  likelihoodRoutine_Data.massVector = &v_mass;
  likelihoodRoutine_Data.timeVector = &v_time;

  QUESO::GenericScalarFunction<QUESO::GslVector,QUESO::GslMatrix>
    likelihoodFunctionObj("cal_like_",paramDomain,likelihoodRoutine,(void *) &likelihoodRoutine_Data,true); // routine computes [ln(function)]
  
  //------------------------------------------------------
  // SIP Step 4 of 6: Define the prior RV
  //------------------------------------------------------
  QUESO::UniformVectorRV<QUESO::GslVector,QUESO::GslMatrix>
    priorRvA("priorA_", paramDomainA);

  QUESO::GslVector meanVector(paramSpaceB.zeroVector() );
  meanVector[0] = input_data("v_l_a_mean",0.001204);
  meanVector[1] = input_data("v_ic_mean",1.0);
  meanVector[2] = input_data("v_l_p_mean",0.001204);
  meanVector[3] = input_data("v_l_k_mean",0.001204);

  QUESO::GslMatrix covMatrix = QUESO::GslMatrix(paramSpaceB.zeroVector());
  covMatrix(0,0) = std::pow(input_data("v_l_a_sigm",0.001009),2);
  covMatrix(1,1) = std::pow(input_data("v_ic_sigm",2.0),2);
  covMatrix(2,2) = std::pow(input_data("v_l_p_sigm",0.001009),2);
  covMatrix(3,3) = std::pow(input_data("v_l_k_sigm",0.001009),2);

  QUESO::GaussianVectorRV<QUESO::GslVector,QUESO::GslMatrix>
    priorRvB("priorB_",paramDomainB,meanVector,covMatrix);

  QUESO::ConcatenatedVectorRV<QUESO::GslVector,QUESO::GslMatrix>
    priorRv("cal_prior_", priorRvA, priorRvB, paramDomain);

  //------------------------------------------------------
  // SIP Step 5 of 6: Instantiate the inverse problem
  //------------------------------------------------------
  // Extra prefix before the default "rv_" prefix
  QUESO::GenericVectorRV<QUESO::GslVector,QUESO::GslMatrix> postRv("cal_post_", paramSpace);

  // No extra prefix before the default "ip_" prefix
  QUESO::StatisticalInverseProblem<QUESO::GslVector,QUESO::GslMatrix> ip("cycle_cal_", NULL, priorRv, likelihoodFunctionObj, postRv);

  //------------------------------------------------------
  // SIP Step 6 of 6: Solve the inverse problem, that is,
  // set the 'pdf' and the 'realizer' of the posterior RV
  //------------------------------------------------------

  ip.solveWithBayesMLSampling();

  //================================================================
  // Statistical forward problem (SFP)
  //================================================================
  cout << endl << endl << "                       Starting Forward Problem Step" << endl << endl;
  // SFP Step 1 of 6: Instantiate the parameter qoi spaces.
  QUESO::VectorSpace<> qoiSpace(env, "qoi_", Size, NULL);

  // SFP Step 2 of 6: Instantiate the parameter domain

  // Not necessary because input RV of the SFP = output RV of SIP.

  // SFP Step 3 of 6: Instantiate the qoi object to be used by QUESO.

  QUESO::GenericVectorFunction<QUESO::GslVector,QUESO::GslMatrix,QUESO::GslVector,QUESO::GslMatrix>
    qoiFunctionObj("qoi_",paramDomain,qoiSpace,qoiRoutine,(void *) &likelihoodRoutine_Data);
  
  // SFP Step 4 of 6: Define the input RV

  // Not necessary because input RV of SFP = output RV of SIP (postRv).

  // SFP Step 5 of 6: Instantiate the forward problem
  QUESO::GenericVectorRV<> qoiRv("qoi_", qoiSpace);

  QUESO::StatisticalForwardProblem<> fp("", NULL, postRv, qoiFunctionObj, qoiRv);

  // SFP Step 6 of 6: Solve the forward problem
  fp.solveWithMonteCarlo(NULL);

}
