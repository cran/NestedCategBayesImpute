#include <Rcpp.h>
using namespace Rcpp;
#include "samplehouseholds.h"

// [[Rcpp::export]]
NumericMatrix samplehouseholds(NumericMatrix phi, NumericMatrix omega, NumericVector pi,
               NumericVector d, List lambda,
               int currrentbatch, int nHouseholds,  int householdsize) {

  int FF = omega.nrow();
  int SS = omega.ncol();
  int p = d.length();
  int n_lambdas = lambda.length();
  int *lambda_columns = new int[n_lambdas];
  double **lambdas = new double*[n_lambdas];
  int maxDDtp = phi.nrow();
  int maxdd = maxDDtp / p;

  //int ncol = householdsize * DIM + 1 + householdsize;
  //output data: zero-based
  //column 0: household unique index
  //column 1: member index within the household (pernum:person number?)
  //column 2 to 2+p-1: individual data
  //column 2+p: 2+p+n_lambdas-2: household level data
  //column 2+p+n_lambdas-1: household group indicator
  //column last hh_size: individual group indicator
  int DIM = 2 + p + n_lambdas - 1; //not output the household size
  int ncol = DIM * householdsize + householdsize  + 1;
  NumericMatrix data(nHouseholds, ncol);

  //copy data from list of matrices to C++
  for (int i = 0; i < n_lambdas; i++) {
    NumericMatrix l = lambda[i];
    lambda_columns[i] = l.ncol();
    lambdas[i] = new double[l.length()];
    std::copy(l.begin(), l.end(), lambdas[i]);
  }
  //printf("in samplehouseholds\n");
  NumericVector rand = runif(nHouseholds * ncol); //at most this many
  sampleHouseholds_imp(data.begin(), rand.begin(), lambdas, lambda_columns, omega.begin(),
                   phi.begin(), pi.begin(),d.begin(),
                   nHouseholds, householdsize, FF, SS,maxdd,p, currrentbatch,n_lambdas);

  //clean up
  delete [] lambda_columns;
  for (int i = 0; i < n_lambdas; i++) {
    delete [] lambdas[i];
  }
  delete [] lambdas;
  //printf("done samplehouseholds\n");
  return data;
}

// [[Rcpp::export]]
NumericMatrix samplehouseholds_HHhead_at_group_level(NumericMatrix phi, NumericMatrix omega, NumericVector pi,
                               NumericVector d, List lambda,
                               int currrentbatch, int nHouseholds,  int householdsize) {

  int FF = omega.nrow();
  int SS = omega.ncol();
  int p = d.length();
  int n_lambdas = lambda.length();
  int *lambda_columns = new int[n_lambdas];
  double **lambdas = new double*[n_lambdas];
  int maxDDtp = phi.nrow();
  int maxdd = maxDDtp / p;

  //int ncol = householdsize * DIM + 1 + householdsize;
  //output data: zero-based
  //column 0: household unique index
  //column 1: member index within the household (pernum:person number?)
  //column 2 to 2+p-1: individual data
  //column 2+p: 2+p+n_lambdas-2: household level data
  //column 2+p+n_lambdas-1: household group indicator
  //column last hh_size: individual group indicator
  int DIM = 2 + p + n_lambdas - 1; //not output the household size
  int ncol = DIM * householdsize + householdsize  + 1;
  NumericMatrix data(nHouseholds, ncol);

  //copy data from list of matrices to C++
  for (int i = 0; i < n_lambdas; i++) {
    NumericMatrix l = lambda[i];
    lambda_columns[i] = l.ncol();
    lambdas[i] = new double[l.length()];
    std::copy(l.begin(), l.end(), lambdas[i]);
  }
  //printf("in samplehouseholds\n");
  NumericVector rand = runif(nHouseholds * ncol); //at most this many
  sampleHouseholds_imp_HHhead_at_group_level(data.begin(), rand.begin(), lambdas, lambda_columns, omega.begin(),
                       phi.begin(), pi.begin(),d.begin(),
                       nHouseholds, householdsize, FF, SS,maxdd,p, currrentbatch,n_lambdas);

  //clean up
  delete [] lambda_columns;
  for (int i = 0; i < n_lambdas; i++) {
    delete [] lambdas[i];
  }
  delete [] lambdas;
  //printf("done samplehouseholds\n");
  return data;
}
