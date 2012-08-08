
#include <R.h>
#include <Rmath.h>
#include "utility.h"

/* Utility function to determine if any element in a vector is negative */
int anyNegative(int n, const int *v)                     // both should be const
{
    int i;
    for (i=0; i<n; i++) 
    {
        if (v[i]<0) return 1;
    }
    return 0;
}


/* Integer vector copy */
void copy(int n, const int *from, int *to) 
{
    int i;
    for (i=0; i<n; i++) to[i]=from[i];
}

void runif_vec(int nLength, const double *adLowerBound, const double *adUpperBound, double *adUnifs)
{
    int i;
    GetRNGstate();
    for (i=0; i<nLength; i++) adUnifs[i] = runif(adLowerBound[i], adUpperBound[i]);
    PutRNGstate();
}

void rpois_vec(int nLength, const double *adMean, int *anPoisson)         
{
    int i;
    GetRNGstate();
    for (i=0; i<nLength; i++) anPoisson[i]=rpois(adMean[i]);
    PutRNGstate(); 
}

void rbeta_vec(int nLength, const double *adParamA, const double *adParamB, double *adBetas)
{
    int i;
    GetRNGstate();
    for (i=0; i<nLength; i++) adBetas[i] = rbeta(adParamA[i], adParamB[i]);
    PutRNGstate();
}

void rgamma_vec(int nLength, const double *adShape, const double *adRate, double *adGammas)
{
    int i;
    GetRNGstate();
    for (i=0; i<nLength; i++) adGammas[i] = rgamma(adShape[i], 1.0/adRate[i]); 
    PutRNGstate();

}



