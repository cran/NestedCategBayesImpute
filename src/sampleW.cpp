//  sampleW.cpp
//  Created by Quanli Wang on 2/20/16.

int samplew(double *p, int n, double d) {
    double dsum;
    int i,k;
    dsum = 0;
    double *myw;
    myw = new double[n];
    for (i = 0; i < n;i++) {
        dsum+=p[i];
    }
    if (dsum <=0 ) {dsum =1;}
    myw[0] = p[0] / dsum;
    for (i = 1; i < n;i++) {
        myw[i] = p[i] / dsum + myw[i-1];
    }

    for(k=0;k < n && d>myw[k];k++)
        ;
    delete [] myw;
    if (k == n) {k = n-1;}
     return k+1;
}

void samplew_multi(double *p, int n, double *d,int howmany) {
    double dsum;
    int i,k;
    dsum = 0;
    double *myw;
    myw = new double[n];
    for (i = 0; i < n;i++) {
        dsum+=p[i];
    }
    if (dsum <=0 ) {dsum =1;}
    myw[0] = p[0] / dsum;
    for (i = 1; i < n;i++) {
        myw[i] = p[i] / dsum + myw[i-1];
    }
    for (int h=0; h < howmany; h++) {
        for(k=0;k < n && d[h]>myw[k];k++)
            ;
        if (k == n) {k = n-1;}
        d[h] = k+1;

    }
    delete [] myw;
}

//this version put results into a different place
void samplew_multi2(double *p, int n, double *d, double* result,int howmany) {
    double dsum;
    int i,k;
    dsum = 0;
    double *myw;
    myw = new double[n];
    for (i = 0; i < n;i++) {
        dsum+=p[i];
    }
    if (dsum <=0 ) {dsum =1;}
    myw[0] = p[0] / dsum;
    for (i = 1; i < n;i++) {
        myw[i] = p[i] / dsum + myw[i-1];
    }
    for (int h=0; h < howmany; h++) {
        for(k=0;k < n && d[h]>myw[k];k++)
            ;
        if (k == n) {k = n-1;}
        result[h] = k+1;
    }
    delete [] myw;
}
