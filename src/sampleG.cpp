#include <Rcpp.h>
using namespace Rcpp;
#include "sampleW.h"

// [[Rcpp::export]]
List sampleG(NumericMatrix phi, NumericMatrix data,
                          NumericMatrix omega, NumericVector pi, NumericVector ni,
                          NumericMatrix HHdata, List lambda) {
  int p = data.nrow();
  int nIndividuals = data.ncol();
  int FF = omega.nrow();
  int SS = omega.ncol();
  int maxdd = phi.nrow() / p;
  int n = ni.length();
  //printf("in sampleG\n");
  std::vector<NumericMatrix> Lambdas;
  for (int i = 0; i < lambda.length(); i++) {
    Lambdas.push_back(lambda[i]);
  }

  //printf("FF = %d, SS = %d, p = %d, maxd = %d, nIndividuals = %d, n=%d\n", FF, SS, p, maxdd, nIndividuals,n);

  NumericVector group(n);
  NumericVector indi(nIndividuals);
  int count = 0;

  NumericVector rand = runif(n);
  //use one-d indexing here to be consistant with Matalb code
  //might need to abandon this if we are going to abondon the Matlab version
  double *Gupdateprob1 = new double[FF];
  int *cum_ni = new int[n];
  cum_ni[0] = 0;
  for (int i = 1; i < n; i++) {
    cum_ni[i] = cum_ni[i-1] + (int)ni[i-1];
  }
  int maxDDtP = maxdd*p;
  for (int h = 0; h < n; h++) {
    for (int k=0; k < FF; k++) {
      try {
        double Gupdateprod = 1.0;
        for (int memberindex=0; memberindex < ni[h]; memberindex++){
          int base = (cum_ni[h]+memberindex)*p; //base for data
          double add = 0.0;
          for (int l=0; l < SS; l++) {
            double phiprod = 1.0;
            int phi_base = (int)(maxDDtP*(k*SS+l));
            for (int j=0; j < p; j++) {
              int u = (int)data[base+j]-1;
              phiprod *= phi[phi_base+j*maxdd+u];
            }
            add += omega[FF*l+k]*phiprod;
          } // closing l++
          Gupdateprod *= add;
        } // closing member++

        for (int hv = 0; hv < Lambdas.size(); hv++) {
          Gupdateprod *= Lambdas[hv][(HHdata[h+hv*n]-1)*FF+k];
        }
        Gupdateprob1[k] = pi[k]*Gupdateprod;
      } catch (...) {
        Gupdateprob1[k] = 0;
      }
    } // closing k++

    group[h] = samplew(Gupdateprob1, FF, rand[h]);
    for (int m=0; m < ni[h];m++) {
      indi[count++] = group[h];
    }
  }
  delete [] cum_ni;
  delete [] Gupdateprob1;
  //printf("done sampleG\n");
  return List::create(Named("G", group), Named("G_Individuals", indi));
}

