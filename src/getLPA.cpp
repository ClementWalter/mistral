#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix getLPA(const NumericVector x, const IntegerVector iv, const int Nppp, const int lmax) {
  const int c = x.length();
  NumericMatrix LPA(Nppp,lmax);
  double fill;
  int cur;
  for (int i = 0; i < Nppp; ++i) {
    fill = NA_REAL;
    cur = iv(i);
    for (int j = c-1; j >= 0; --j) {
      if (iv(j) == cur) fill = x(j);
      if(j<lmax) LPA(i,j) = fill;
    }
  }
  return LPA;
}
