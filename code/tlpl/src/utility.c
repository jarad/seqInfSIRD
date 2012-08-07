
#include <R.h>
#include <Rmath.h>
#include "utility.h"

/* Utility function to determine if any element in a vector is negative */
int anyNegative(int n, int *v)                     // both should be const
{
    int i;
    for (i=0; i<n; i++) 
    {
        if (v[i]<0) return 1;
    }
    return 0;
}


/* Integer vector copy */
void copy(int n, int *from, int *to) 
{
    int i;
    for (i=0; i<n; i++) to[i]=from[i];
}


void rpois_vec(const int *nLength, const double *adMean, int *anPoisson)         
{
    int i;
    GetRNGstate();
    for (i=0; i<*nLength; i++) anPoisson[i]=rpois(adMean[i]);
    PutRNGstate(); 
}

void rbeta_vec(const int *nLength, const double *adParamA, const double *adParamB, double *adBetas)
{
    int i;
    GetRNGstate();
    for (i=0; i<*nLength; i++) adBetas[i] = rbeta(adParamA[i], adParamB[i]);
    PutRNGstate();
}

void rgamma_vec(const int *nLength, const double *adShape, const double *adRate, double *adGammas)
{
    int i;
    GetRNGstate();
    for (i=0; i<*nLength; i++) adGammas[i] = rgamma(adShape[i], 1.0/adRate[i]); 
    PutRNGstate();

}



