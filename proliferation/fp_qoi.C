#include <fp_qoi.h>
#include <ip_likelihood.h>
#include <general_libraries.h>
#include <cmath>
#include "main_model.h"

void qoiRoutine(
  const QUESO::GslVector&                    paramValues,
  const QUESO::GslVector*                    paramDirection,
  const void*                                functionDataPtr,
        QUESO::GslVector&                    qoiValues,
        QUESO::DistArray<QUESO::GslVector*>* gradVectors,
        QUESO::DistArray<QUESO::GslMatrix*>* hessianMatrices,
        QUESO::DistArray<QUESO::GslVector*>* hessianEffects)
{
  // Logic just to avoid warnings from INTEL compiler
  const QUESO::GslVector* aux1 = paramDirection;
  if (aux1) {};
  QUESO::DistArray<QUESO::GslVector*>* aux2 = gradVectors;
  if (aux2) {};
  aux2 = hessianEffects;
  QUESO::DistArray<QUESO::GslMatrix*>* aux3 = hessianMatrices;
  if (aux3) {};
  
  unsigned int N_par=5;

  UQ_FATAL_TEST_MACRO(paramValues.sizeGlobal() != N_par,QUESO::UQ_UNAVAILABLE_RANK,"likelihoodRoutine()","paramValues vector does not have the correct size");

  double stdDevs = paramValues[2];
  double icon = paramValues[4];

  std::vector<double> Parameters(3,0);
  Parameters[0] = paramValues[0];
  Parameters[1] = paramValues[1];
  Parameters[2] = paramValues[3];
  
  const vector<double>* timeVector = ((likelihoodRoutine_DataType *) functionDataPtr)->timeVector;
  
  const unsigned int Size = timeVector->size();
  vector<double> out_mod;

  main_code(out_mod,Parameters,icon,(*timeVector));

  const QUESO::BaseEnvironment& env = paramValues.env();
  if (env.subRank() == 0) {
    for (unsigned int i = 0; i < Size; ++i) {
      qoiValues[i] = out_mod[i];
    }
  }
  else {
    // Do nothing;
  }

}
