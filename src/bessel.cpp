#include "rcpplist.h"

double besselJ(double x, double nu)
{
  return boost::math::cyl_bessel_j(nu, x);
}

double besselJ_modified_C(const double x, const double nu)
{
  if(x >= 0.0) {

    return besselJ(x, nu);

  } else {

    return pow(-1.0, nu) * besselJ(-x, nu);
  }
}
