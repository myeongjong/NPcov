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

double est_2d_epan_refl_C(const double r, const Rcpp::NumericVector data, const double h)
{
  double value = 0.0;
  for(unsigned int j = 0; j < data.length(); j++) {

    value += (1-pow(data[j]/h, 2.0)) * (lmda1_C(r*(data[j]+h)) - lmda1_C(r*(data[j]-h))) - ((pow(r,2.0)*pow(data[j]+h, 2.0) * besselJ_modified_C(r*(data[j]+h), 1.0) - lmda0_C(r*(data[j]+h))) - (pow(r, 2.0)*pow(data[j]-h, 2.0) * besselJ_modified_C(r*(data[j]-h), 1.0) - lmda0_C(r*(data[j]-h)))) / pow(r * h, 2.0) + (2 * data[j] / r / pow(h, 2.0)) * (r*(data[j]+h) * besselJ_modified_C(r*(data[j]+h), 1.0) - r*(data[j]-h) * besselJ_modified_C(r*(data[j]-h), 1.0));
  }

  return value * 3.0 / 4.0 / (double)data.length() / h / r;
}

double est_2d_gaussian_refl_C(const double r, const Rcpp::NumericVector data, const double h)
{
  double value = 0.0;
  for(unsigned int j = 0; j < data.length(); j++) {

    value += besselJ_generalized_C(r, h, data[j]);
  }

  return value * exp(-0.25 * pow(h, 2.0) * pow(r, 2.0)) / M_PI / (double)data.length();
}

double est_2d_uniform_refl_C(const double r, const Rcpp::NumericVector data, const double h)
{
  double value = 0.0;
  for(unsigned int j = 0; j < data.length(); j++) {

    value += lmda1_C(r*(data[j]+h)) - lmda1_C(r*(data[j]-h));
  }

  return value / 2.0 / (double)data.length() / h / r;
}

double est_infd_epan_refl_C(const double r, const Rcpp::NumericVector data, const double h)
{
  double value = 0.0;
  for(unsigned int j = 0; j < data.length(); j++) {

    value += (1-pow(data[j]/h, 2.0)) * (M_SQRT_PI/r) * (R::pnorm(M_SQRT2 * r * (data[j]+h), 0.0, 1.0, true, false) - R::pnorm(M_SQRT2 * r * (data[j]-h), 0.0, 1.0, true, false)) + (data[j] / pow(h, 2.0) / pow(r, 2.0)) * (exp(-pow(r * (data[j]-h), 2.0)) - exp(-pow(r * (data[j]+h), 2.0))) - (0.5 / pow(r, 2.0) / pow(h, 2.0)) * (((data[j]-h)*exp(-pow(r * (data[j]-h), 2.0)) - (data[j]+h)*exp(-pow(r * (data[j]+h), 2.0))) + (M_SQRT_PI/r) * (R::pnorm(M_SQRT2 * r * (data[j]+h), 0.0, 1.0, true, false) - R::pnorm(M_SQRT2 * r * (data[j]-h), 0.0, 1.0, true, false)));
  }

  return value * 3.0 / 4.0 / (double)data.length() / h;
}

double est_infd_gaussian_refl_C(const double r, const Rcpp::NumericVector data, const double h)
{
  double value = 0.0, A = pow(r, 2.0) + 1 / (2 * pow(h, 2.0));
  for(unsigned int j = 0; j < data.length(); j++) {

    value += exp(pow(data[j]/pow(h, 2.0), 2.0) / 4.0 / A - pow(data[j], 2.0) / 2.0 / pow(h, 2.0));
  }

  return value / ((double)data.length() * h) / M_SQRT2 / pow(A, 0.5);
}

double est_infd_uniform_refl_C(const double r, const Rcpp::NumericVector data, const double h)
{
  double value = 0.0;
  for(unsigned int j = 0; j < data.length(); j++) {

    value += R::pnorm(M_SQRT2 * r *(data[j]+h), 0.0, 1.0, true, false) - R::pnorm(M_SQRT2 * r *(data[j]-h), 0.0, 1.0, true, false);
  }

  return value * M_SQRT_PI / 2.0 / (double)data.length() / h / r;
}

double estimator_2d_reflection(const double r, const Rcpp::NumericVector data, const double h, const std::string kernel)
{
  if(kernel == "epan") {

    return est_2d_epan_refl_C(r, data, h);

  } else if(kernel == "gaussian") {

    return est_2d_gaussian_refl_C(r, data, h);

  } else if(kernel == "uniform"){

    return est_2d_uniform_refl_C(r, data, h);

  } else {

    return NA_REAL;
  }
}

double estimator_infd_reflection(const double r, const Rcpp::NumericVector data, const double h, const std::string kernel)
{
  if(kernel == "epan") {

    return est_infd_epan_refl_C(r, data, h);

  } else if(kernel == "gaussian") {

    return est_infd_gaussian_refl_C(r, data, h);

  } else if(kernel == "uniform"){

    return est_infd_uniform_refl_C(r, data, h);

  } else {

    return NA_REAL;
  }
}

//' @name ksreflect_C
//'
//' @title Kernel-type estimator based on the reflection method (C++ version)
//'
//' @param r Distance
//' @param data Pseudo data
//' @param h Bandwidth
//' @param kernel Type of kernel: epan, gaussian, or uniform
//' @param shape Shape of estimator: general or monotone
//'
//' @return Numeric value of the estimator
// [[Rcpp::export]]
double ksreflect_C(const double r, const Rcpp::NumericVector data, const double h, const std::string kernel, const std::string shape)
{
  if(shape == "general") {

    return estimator_2d_reflection(r, data, h, kernel);

  } else if(shape == "monotone") {

    return estimator_infd_reflection(r, data, h, kernel);

  } else {

    return NA_REAL;
  }
}

//' @name eval_sse_C
//'
//' @title Computing SSE of an estimator obtained by the IDEA algorithm (C++ version)
//'
//' @param storage Data storage
//' @param x Input (e.g., distance)
//' @param y Output (e.g., correlation value)
//' @param h Bandwidth
//' @param kernel Type of kernel: epan, gaussian, or uniform
//' @param shape Shape of estimator: general or monotone
//'
//' @return SSE
// [[Rcpp::export]]
Rcpp::NumericVector eval_sse_C(const Rcpp::NumericMatrix storage, const Rcpp::NumericVector x, const Rcpp::NumericVector y, const double h, const std::string kernel, const std::string shape)
{
  Rcpp::NumericVector tempvec (storage.nrow());
  Rcpp::NumericVector obj_value (storage.ncol());
  for(unsigned int j = 0; j < storage.ncol(); j++) {

    tempvec = storage.column(j);
    for(unsigned int k = 0; k < x.length(); k++) {

      obj_value[j] += pow(y[k] - ksreflect_C(x[k], tempvec, h, kernel, shape), 2.0);
    }
  }

  return obj_value;
}

//' @name evalEach_sse_C
//'
//' @title Computing SSE of an estimator obtained by the IDEA algorithm for each pseudo data (C++ version)
//'
//' @param data Pseudo data
//' @param x Input (e.g., distance)
//' @param y Output (e.g., correlation value)
//' @param h Bandwidth
//' @param kernel Type of kernel: epan, gaussian, or uniform
//' @param shape Shape of estimator: general or monotone
//'
//' @return SSE
// [[Rcpp::export]]
double evalEach_sse_C(const Rcpp::NumericVector data, const Rcpp::NumericVector x, const Rcpp::NumericVector y, const double h, const std::string kernel, const std::string shape)
{
  double obj_value = 0.0;
  for(unsigned int k = 0; k < x.length(); k++) {

    obj_value += pow(y[k] - ksreflect_C(x[k], data, h, kernel, shape), 2.0);
  }

  return obj_value;
}
