#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include "utility.h"

/*********************** Utility functions - to be moved *******************************/


int double_comp(const void *X, const void *Y) 
{
    double x = *((double *)X);
    double y = *((double *)Y);
    if (x>y) 
    {
        return 1;
    }
    else 
    {
        if (x<y)
        {
            return -1;
        }
        else
        {
            return 0;
        }
    }
}

int is_sorted(const int n, const double *v) 
{
    int i; 
    for (i=1; i<n; i++)
    {
        if (v[i]<v[i-1]) return 0;
    }
    return 1;
}

void cumulative_sum(const int n, double *v) 
{
    int i;
    for (i=1; i<n; i++) v[i] += v[i-1];
}



/************************ Resampling functions *****************************************/

void inverse_cdf_weights(const int *nWeights, double *adWeights, const int *nUniforms, const double *adUniforms,
                         int *anIndices)
{
    if (!is_sorted(*nUniforms, adUniforms))    
        qsort((void *)adUniforms, *nUniforms, sizeof(double), double_comp);

    cumulative_sum(*nWeights, adWeights);

    int i, j=0, found;
    for (i=0; i<*nUniforms; i++) 
    {
        found=0;
        while (!found) 
        {
            if (adUniforms[i] > adWeights[j])
            {
               j++;
            }
            else 
            {
                found=1;
            }
        }
        anIndices[i] = j;
    }    
}






/* Renormalizes the vector adWeights to sum to 1. If *nLog \ne 0, the weights are assumed to be log weights */
void renormalize(const int *nWeights, const int *nLog, double *adWeights) 
{
    int i;

    if (*nLog) {
        double max=adWeights[0];
        for (i=1; i<*nWeights; i++) if (adWeights[i]>max) max = adWeights[i]; 
        for (i=0; i<*nWeights; i++) adWeights[i] = exp(adWeights[i]-max);
    }

    double sum=0;
    for (i=0; i<*nWeights; i++) sum += adWeights[i];
    for (i=0; i<*nWeights; i++) adWeights[i] /= sum;
}


void multinomial_resample(const int *nWeights, double *adWeights, const int *nIndices, int *anIndices) 
{
    int i, j;
    double cusum[*nWeights], adUniforms[*nIndices];
    cumulative_sum(*nWeights, adWeights);

    GetRNGstate();
    for (i=0; i<*nIndices; i++) adUniforms[i] = runif(0,1);
    PutRNGstate();

    inverse_cdf_weights(nWeights, adWeights, nIndices, adUniforms, anIndices);
}


void stratified_resample(const int *nWeights, double *adWeights, const int *nIndices, int *anIndices)
{
    int i;
    double adLowerBound[*nIndices], adUpperBound[*nIndices];
    for (i=0;i<*nIndices;i++)
    {
        adLowerBound[i] = (float)  i    / *nIndices;
        adUpperBound[i] = (float) (i+1) / *nIndices;
    }
    
    double adUniforms[*nIndices];
    runif_vec(nIndices, adLowerBound, adUpperBound, adUniforms);

    Rprintf("LowerBound: "); for (i=0;i<*nIndices;i++) Rprintf("%1.6f ", adLowerBound[i]); Rprintf("\n");
    Rprintf("UpperBound: "); for (i=0;i<*nIndices;i++) Rprintf("%1.6f ", adUpperBound[i]); Rprintf("\n");  
    Rprintf("Uniforms: "); for (i=0;i<*nIndices;i++) Rprintf("%1.6f ", adUniforms[i]); Rprintf("\n"); 
    inverse_cdf_weights(nWeights, adWeights, nIndices, adUniforms, anIndices);
}



