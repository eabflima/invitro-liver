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

  //unsigned int N_par=4;
  
  //UQ_FATAL_TEST_MACRO(paramValues.sizeGlobal() != N_par,QUESO::UQ_UNAVAILABLE_RANK,"likelihoodRoutine()","paramValues vector does not have the correct size");

  std::vector<double> Parameters(6,0);
  
  Parameters[0] = paramValues[0];
  Parameters[1] = paramValues[1];
  Parameters[2] = paramValues[2];
  Parameters[3] = 0.0;
  Parameters[4] = 1.0;
  Parameters[5] = 0.0;
  
  double stdDevs = paramValues[3];

  vector<double> out_rd;
  
  int* argc   = ((likelihoodRoutine_DataType *) functionDataPtr)->argc;
  char** argv = ((likelihoodRoutine_DataType *) functionDataPtr)->argv;
  const QUESO::BaseEnvironment& env = paramValues.env();
  MPI_Comm libMeshComm = env.subComm().Comm();

  main_code(libMeshComm, argc[0], argv, out_rd, Parameters);

  double result = 0.;
  const unsigned int Size = out_rd.size();

  if (env.subRank() == 0) {
    for (unsigned int i = 0; i < Size; i++) {
      double variance = pow(stdDevs,2);    
      double ratio = (1.0-out_rd[i]);
      result += (std::pow(ratio,2)/(2.0*variance))+log(stdDevs)+0.5*log(M_PI)+0.5*log(2.0);//Gaussian
    }
  }
  else {
    // Do nothing;
  }
  
  return -1.0*result;
  
}
