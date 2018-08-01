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

#include "libmesh/libmesh.h"
//#include "libmesh/getpot.h"

using namespace libMesh;
using namespace std;

void computeIP(const QUESO::FullEnvironment& env, int argc, char** argv)
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
  //GetPot input_data("inp_dat.in");

  unsigned int N_par=4;

  QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix> paramSpace(env, "param_", N_par, NULL);

  //------------------------------------------------------
  // SIP Step 2 of 6: Instantiate the parameter domain
  //------------------------------------------------------
  QUESO::GslVector paramMinValues(paramSpace.zeroVector());
  QUESO::GslVector paramMaxValues(paramSpace.zeroVector());
  
  paramMinValues[0] = 1e-4;//input_data("v_m_v_min",0.2);
  paramMaxValues[0] = 1000.0;//input_data("v_m_v_max",0.10);
  paramMinValues[1] = 0.15;//input_data("v_E_t_min",1000.0);
  paramMaxValues[1] = 1.15;//input_data("v_E_t_max",200000.0);
  paramMinValues[2] = 3.0;//input_data("v_eps_min",10.0);
  paramMaxValues[2] = 1000.0;//input_data("v_eps_max",20000.0);
  paramMinValues[3] = 1e-4;//stdev
  paramMaxValues[3] = 0.2;//stdev

  QUESO::BoxSubset<QUESO::GslVector,QUESO::GslMatrix> paramDomain("param_",paramSpace,paramMinValues,paramMaxValues);

  //------------------------------------------------------
  // SIP Step 3 of 6: Instantiate the likelihood function
  // object to be used by QUESO.
  //------------------------------------------------------

  likelihoodRoutine_DataType likelihoodRoutine_Data;
  int *argp = new int[1];
  argp[0] = argc;
  likelihoodRoutine_Data.argc = argp;
  likelihoodRoutine_Data.argv = argv;

  QUESO::GenericScalarFunction<QUESO::GslVector,QUESO::GslMatrix>
    likelihoodFunctionObj("cal_like_",paramDomain,likelihoodRoutine,(void *) &likelihoodRoutine_Data,true); // routine computes [ln(function)]
  
  //------------------------------------------------------
  // SIP Step 4 of 6: Define the prior RV
  //------------------------------------------------------
  QUESO::UniformVectorRV<QUESO::GslVector,QUESO::GslMatrix> priorRv("cal_prior_", paramDomain);
    
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
}
