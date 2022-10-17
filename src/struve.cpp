/*
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 Copyright (C) 2015 -- 2021 Martin Schlather

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 3
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#include "rcpplist.h"
using namespace Rcpp;

#define RF_NA NA_REAL

#define EXP std::exp
#define LOG std::log
#define POW(X, Y) R_pow((double) X, (double) Y) // OK; keine Klammern um X!
#define FABS(X) std::fabs((double) X) // OK; keine Klammern um X!
#define LENGTH length // to avoid the unvoluntiered use of LENGTH defined by R

double struve_intern(double x, double nu, double factor_Sign, bool expscaled)
{
  if ((x == 0.0) && (nu>-1.0)) return 0.0;
  if (x <= 0.0) return RF_NA; // not programmed yet
  double exp_dummy,
  dummy = 0.0,
  logx = 2.0 * LOG(0.5 * x),
  x1 = 1.5,
  x2 = nu + 1.5,
  value = 1.0,
  fsign = factor_Sign,
  epsilon=1e-20;

  do {
    dummy += logx - LOG(x1) - LOG(FABS(x2));
    exp_dummy = EXP(dummy);
    value += (1 - 2 * (x2 < 0))  * fsign * exp_dummy;
    //  printf("%10g %10g %10g %10g\n", value, fsign, x1, x2);
    x1 += 1.0;
    x2 += 1.0;
    fsign = factor_Sign * fsign;
  } while (exp_dummy > FABS(value) * epsilon);

  x1 = 1.5;
  x2 = nu + 1.5;
  if (x2 > 0.0) {
    dummy = (nu + 1.0) * 0.5 * logx - R::lgammafn(x1) - R::lgammafn(x2);
    if (expscaled) dummy -= x;
    value *= EXP(dummy);
  } else {
    //if ( (double) ((int) (x1-0.5)) != x1-0.5 ) return RF_NA;
    value *= POW(0.5 * x, nu + 1.0) / (R::gammafn(x1) * R::gammafn(x2));
    if (expscaled) value *= EXP(-x);
  }

  return value;
}

double struveH(const double x, const double nu)
{
  return struve_intern(x, nu, -1.0, false);
}

double struveL(const double x, const double nu, const bool expScaled)
{
  return struve_intern(x, nu, 1.0, expScaled);
}

double struveH_modified_C(const double x, const double nu)
{
  if(x >= 0.0) {

    return struveH(x, nu);

  } else {

    return pow(-1.0, nu + 1) * struveH(x, nu);
  }
}

