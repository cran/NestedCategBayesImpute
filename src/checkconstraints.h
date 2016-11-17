//  sampleHouseholds.h
#define HEAD 1
#define SPOUSE 2
#define CHILD 3
#define CHILDINLAW 4
#define PARENT 5
#define PARENTINLAW 6
#define SIBLING 7
#define SIBLINGINLAW 8
#define GRANDCHILD 9

#define COL 3

int checkconstraints_imp(double *data, double *isPossible,int hh_size, int DIM, int nHouseholds);
int checkconstraints_imp_HHhead_at_group_level(double *data, double *isPossible,int hh_size, int DIM, int nHouseholds);
