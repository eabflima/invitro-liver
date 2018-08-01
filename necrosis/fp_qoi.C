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

  GetPot input_data("inp_dat.in");
  const double nut_mean        = input_data("v_nut_val",1.0);
  
  unsigned int N_par=6;

  UQ_FATAL_TEST_MACRO(paramValues.sizeGlobal() != N_par,QUESO::UQ_UNAVAILABLE_RANK,"likelihoodRoutine()","paramValues vector does not have the correct size");

  double stdDevs = paramValues[1];
  double icon = paramValues[3];

  std::vector<double> Parameters(4,0);
  Parameters[0] = paramValues[4]; //prol
  Parameters[1] = paramValues[5]; //supp
  Parameters[2] = paramValues[2]; //apop
  Parameters[3] = paramValues[0]; //nec

  const vector<double>* timeVector = ((likelihoodRoutine_DataType *) functionDataPtr)->timeVector;
  
  const unsigned int Size = timeVector->size();
  vector<double> out_mod;

  main_code(out_mod,Parameters,icon,nut_mean,(*timeVector));

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
