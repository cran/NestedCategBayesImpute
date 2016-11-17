#include <Rcpp.h>
using namespace Rcpp;
#include "sampleW.h"

// [[Rcpp::export]]
NumericVector sampleM(NumericMatrix phi, NumericMatrix data,
               NumericMatrix omega, NumericVector G, NumericVector serial) {

  //printf("in sampleM\n");
  int p = data.nrow();
  int nIndividuals = data.ncol();
  int FF = omega.nrow();
  int SS = omega.ncol();
  int maxDDtp = phi.nrow();
  int maxdd = maxDDtp / p;
  //printf("maxdd = %d\n", maxdd);
  NumericVector indi(nIndividuals);
  NumericVector rand = runif(nIndividuals);

  double *Gupdateprob2= new double[SS];
  for (int m = 0; m < nIndividuals; m++) {

    int G_asg = G[int(serial[m])-1];
    int base = m*p;
    for (int l = 0; l < SS; l++) {
      try {
        double phiprod = 1.0;
        for (int j = 0; j < p; j++) {
          int u = (int)data[base+j]-1;
          phiprod *= phi[maxDDtp*((G_asg-1)*SS+l)+j*maxdd+u];
        }
        Gupdateprob2[l] = omega[FF*l+G_asg-1]*phiprod;
      } catch(...) {
        Gupdateprob2[l] = 0;
      }
    }
    indi[m] = samplew(Gupdateprob2, SS, rand[m]);
    //printf("indi = (%d,%d)\n",m, (int)indi[m]);
  }
  delete [] Gupdateprob2;
  //printf("done sampleM\n");
  return indi;
}

