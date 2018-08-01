#include <ip_likelihood.h>
#include <general_libraries.h>
#include "main_model.h"

double likelihoodRoutine(
  const QUESO::GslVector& paramValues,
  const QUESO::GslVector* paramDirection,
  const void*             functionDataPtr,
  QUESO::GslVector*       gradVector,
  QUESO::GslMatrix*       hessianMatrix,
  QUESO::GslVector*       hessianEffect)
{
  // Logic just to avoid warnings from INTEL compiler
  const QUESO::GslVector* aux1 = paramDirection;
  if (aux1) {};
  aux1 = gradVector;
  aux1 = hessianEffect;
  QUESO::GslMatrix* aux2 = hessianMatrix;
  if (aux2) {};

  GetPot input_data("inp_dat.in");
  const unsigned int n_samples = input_data("n_samples",4);

  unsigned int N_par=5;
  
  UQ_FATAL_TEST_MACRO(paramValues.sizeGlobal() != N_par,QUESO::UQ_UNAVAILABLE_RANK,"likelihoodRoutine()","paramValues vector does not have the correct size");

  double stdDevs = paramValues[2];
  double icon = paramValues[4];

  std::vector<double> Parameters(3,0);
  Parameters[0] = paramValues[0];
  Parameters[1] = paramValues[1];
  Parameters[2] = paramValues[3];

  const vector<double>* massVector = ((likelihoodRoutine_DataType *) functionDataPtr)->massVector;
  const vector<double>* timeVector = ((likelihoodRoutine_DataType *) functionDataPtr)->timeVector;

  const unsigned int Size = timeVector->size();
  vector<double> out_mod;

  main_code(out_mod,Parameters,icon,(*timeVector));

  double result = 0.;
  const QUESO::BaseEnvironment& env = paramValues.env();
  if (env.subRank() == 0) {
    for (unsigned int i = 0; i < Size; i++) {
      double variance = pow(stdDevs,2);    
      double out_rd = out_mod[i];
      for (unsigned int j = 0; j < n_samples; j++) {
	/*
	//////====== New Log-Gaussian
	double u = log(out_rd);
	double x = ((*massVector)[i*n_samples+j]);
	double s = variance;
	result += 0.5*pow(log(x),2)/s
	  -u*log(x)/s
	  +log(x)
	  +log(sqrt(s))
	  +0.5*pow(u,2)/s
	  +0.5*log(M_PI)
	  +0.5*log(2.0);//Log-Gaussian
	*/
	double ratio = (out_rd - (*massVector)[i*n_samples+j]);
	result += (std::pow(ratio,2)/(2.0*variance))+log(stdDevs)+0.5*log(M_PI)+0.5*log(2.0);//Gaussian
	/*
	double x = out_rd;
	double m = (*massVector)[i*n_samples+j];
	//double v = variance;
	double s = log(1.0+(stdDevs/pow(m,2)));
	double u = log(m/sqrt(1.0+(stdDevs/pow(m,2))));
	result += 0.5*pow(log(x),2)/s
	  -u*log(x)/s
	  +log(x)
	  +log(sqrt(s))
	  +0.5*pow(u,2)/s
	  +0.5*log(M_PI)
	  +0.5*log(2.0);//Log-Gaussian
	*/
      }
    }
  }
  else {
    // Do nothing;
  }

  return -1.0*result;

}
