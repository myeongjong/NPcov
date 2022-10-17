#include "rcpplist.h"

double lmda0_C(const double x)
{
  double output = (M_PI * x / 2.0) * (besselJ_modified_C(x, 1.0) * struveH_modified_C(x, 0.0) - besselJ_modified_C(x, 0.0) * struveH_modified_C(x, 1.0));

  return output;
}

double lmda1_C(const double x)
{
  double output = x * besselJ_modified_C(x, 0.0) + lmda0_C(x);

  return output;
}
