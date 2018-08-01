#ifndef EX_LIKELIHOOD_H
#define EX_LIKELIHOOD_H

#include <queso/GslMatrix.h>

struct
likelihoodRoutine_DataType
{
  const std::vector<double>* massVector;
  const std::vector<double>* timeVector;
};

double likelihoodRoutine(
  const QUESO::GslVector& paramValues,
  const QUESO::GslVector* paramDirection,
  const void*             functionDataPtr,
  QUESO::GslVector*       gradVector,
  QUESO::GslMatrix*       hessianMatrix,
  QUESO::GslVector*       hessianEffect);

#endif
