//  sampleHouseholds.h
//#define DIM 8
void sampleHouseholds_imp(double* data, double* rand,  double** lambda, int* lambda_columns, double* omega, double* phi,
                      double *pi, double* d,int nHouseholds, int householdsize, int FF,int SS,
                      int maxdd, int p, int currrentbatch, int n_lambdas);
void sampleHouseholds_imp_HHhead_at_group_level(double* data, double* rand,  double** lambda, int* lambda_columns, double* omega, double* phi,
                          double *pi, double* d,int nHouseholds, int householdsize, int FF,int SS,
                          int maxdd, int p, int currrentbatch, int n_lambdas);
