#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix groupcount(NumericVector g1, NumericVector g2, int n1, int n2) {

  NumericMatrix counts(n1, n2);
  for (int i = 0; i < g1.length(); i++) {
    counts[((int)g1[i] -1) + ((int)g2[i]-1) * n1 ]++;
  }
  return counts;
}

// [[Rcpp::export]]
NumericVector groupcount1D(NumericVector g, int n) {
  NumericVector counts(n);
  for (int i = 0; i < g.length(); i++) {
    counts[((int)g[i]-1)]++;
  }
  return counts;
}
