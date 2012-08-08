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

void rep2id(const int *n, int *rep, int *sum, int *id)
{
    int i, j=0;

    i=0;
    while (i<*sum) 
    {
        if (rep[j]>0) // If this particle is resampled (again)
        {
            rep[j]--; 
            id[i] = j;
            i++;
        } 
        else          // If all resamples of this particle are exhausted
        {
            j++;
        }
    }
}


void inverse_cdf_weights(const int *nWeights, 
                         double *adWeights, 
                         const int *nUniforms, 
                         const double *adUniforms,
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





/************************ Resampling functions *****************************************/








/* Renormalizes the vector adWeights to sum to 1. If *nLog \ne 0, the weights are assumed to be log weights */
void renormalize(int nWeights, int nLog, double *adWeights) 
{
    int i;

    if (nLog) {
        double max=adWeights[0];
        for (i=1; i<nWeights; i++) if (adWeights[i]>max) max = adWeights[i]; 
        for (i=0; i<nWeights; i++) adWeights[i] = exp(adWeights[i]-max);
    }

    double sum=0;
    for (i=0; i<nWeights; i++) sum += adWeights[i];
    for (i=0; i<nWeights; i++) adWeights[i] /= sum;
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

    inverse_cdf_weights(nWeights, adWeights, nIndices, adUniforms, anIndices);
}


void systematic_resample(const int *nWeights, double *adWeights, const int *nIndices, int *anIndices)
{
    int i;
    double adUniforms[*nIndices];
    GetRNGstate();
    adUniforms[0] = runif(0, (float) 1/ *nIndices);
    PutRNGstate();
    for (i=1; i<*nIndices; i++) adUniforms[i] =  adUniforms[i-1]+ (float) 1 / *nIndices;

    inverse_cdf_weights(nWeights, adWeights, nIndices, adUniforms, anIndices);
}

void residual_resample(const int *nWeights, double *adWeights, int *nIndices, int *anIndices,
                       const int *nResidualResampleFunction)
{
    int i, anDeterministicReps[*nWeights], nDeterministicReps=0;
    double adExpectedSamples[*nWeights];
    for (i=0; i<*nWeights; i++) 
    {
        adExpectedSamples[i]    = adWeights[i]* *nIndices;
        anDeterministicReps[i]  = adExpectedSamples[i];
        adWeights[i]            = adExpectedSamples[i]-anDeterministicReps[i];
        nDeterministicReps     += anDeterministicReps[i];
    }
    if (nDeterministicReps > *nIndices) 
        REprintf("C: residual_resample: too many deterministic reps\n");
   
    rep2id(nWeights, anDeterministicReps, &nDeterministicReps, anIndices);
    *nIndices -= nDeterministicReps;

    Rprintf("Remaining samples: %d\n", *nIndices);
    Rprintf("Deterministic reps: %d\n", nDeterministicReps);

    renormalize(*nWeights, 0, adWeights);
    switch (*nResidualResampleFunction) 
    {
        case 1:
            stratified_resample( nWeights, adWeights, nIndices, &anIndices[nDeterministicReps]);
            break;
        case 2:
            multinomial_resample(nWeights, adWeights, nIndices, &anIndices[nDeterministicReps]);
            break;
        case 3:
            systematic_resample( nWeights, adWeights, nIndices, &anIndices[nDeterministicReps]);
            break;
        default:
            REprintf("C: residual_resample: no match for residual resampling function\n");
    }
       
}




