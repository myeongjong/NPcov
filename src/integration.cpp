// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

#include "rcpplist.h"
using namespace Numer;

class integrand: public Func
{
private:
  double r;
  double h;
  double t;
public:
  integrand(double r_, double h_, double t_) : r(r_), h(h_), t(t_) {}

  double operator()(const double& x) const
  {
    return exp(0.25 * pow(h, 2.0) * pow(r, 2.0) * cos(2 * x)) * cos(r * t * sin(x));
  }
};

Rcpp::List integrate_list(const double& r, const double& h, const double& t)
{
  const double lower = 0.0, upper = M_PI;

  integrand f(r, h, t);
  double err_est;
  int err_code;
  const double res = integrate(f, lower, upper, err_est, err_code);

  return Rcpp::List::create(
    Rcpp::Named("approximate") = res,
    Rcpp::Named("error_estimate") = err_est,
    Rcpp::Named("error_code") = err_code
  );
}

double besselJ_generalized_C(const double& r, const double& h, const double& t)
{
  const double lower = 0.0, upper = M_PI;

  integrand f(r, h, t);
  double err_est;
  int err_code;
  const double res = integrate(f, lower, upper, err_est, err_code);

  return res;
}
