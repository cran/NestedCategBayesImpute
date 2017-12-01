//  checkconstraints.h
#define HEAD 1
#define SPOUSE 2
#define BIOLOGICALCHILD 3
#define ADOPTEDCHILD 4
#define STEPCHILD 5
#define SIBLING 6
#define PARENT 7
#define GRANDCHILD 8
#define PARENTINLAW 9
#define CHILDINLAW 10

#define COL 3

#define GENDER 0
#define AGE 3
#define RELATE 4

int checkconstraints_imp(double *data, double *isPossible,int hh_size, int DIM, int nHouseholds);
int checkconstraints_imp_HHhead_at_group_level(double *data, double *isPossible,int hh_size, int DIM, int nHouseholds);
