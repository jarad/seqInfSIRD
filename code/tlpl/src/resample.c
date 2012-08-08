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

int is_sorted(int n, const double *v) 
{
    int i; 
    for (i=1; i<n; i++)
    {
        if (v[i]<v[i-1]) return 0;
    }
    return 1;
}

void cumulative_sum(int n, double *v) 
{
    int i;
    for (i=1; i<n; i++) v[i] += v[i-1];
}

void rep2id(int *rep, int sum, int *id)
{
    // This implementation seems poor. 
    // No error checking to assure we stay within the bounds of rep and id
    int i, j=0;

    i=0;
    while (i<sum) 
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


void inverse_cdf_weights(int nW, 
                         double *adWeights, 
                         int nU, 
                         const double *adUniforms,
                         int *anIndices)
{
    if (!is_sorted(nU, adUniforms))    
        qsort((void *)adUniforms, nU, sizeof(double), double_comp);

    cumulative_sum(nW, adWeights);

    int i, j=0, found;
    for (i=0; i<nU; i++) 
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




/* Renormalizes the vector w to sum to 1. 
 * If log \ne 0, the weights are assumed to be log weights */
void renormalize(int n, int log, double *w) 
{
    int i;

    if (log) {
        double max=w[0];
        for (i=1; i<n; i++) if (w[i]>max) max = w[i]; 
        for (i=0; i<n; i++) w[i] = exp(w[i]-max);
    }

    double sum=0;
    for (i=0; i<n; i++) sum += w[i];
    for (i=0; i<n; i++) w[i] /= sum;
}


void multinomial_resample(int nW, double *adWeights, int nI, int *anIndices) 
{
    int i, j;
    double cusum[nW], adUniforms[nI];

    GetRNGstate();
    for (i=0; i<nI; i++) adUniforms[i] = runif(0,1);
    PutRNGstate();

    inverse_cdf_weights(nW, adWeights, nI, adUniforms, anIndices);
}


void stratified_resample(int nW, double *adWeights, int nI, int *anIndices)
{
    int i;
    double adLowerBound[nI], adUpperBound[nI];
    for (i=0;i<nI;i++)
    {
        adLowerBound[i] = (float)  i    / nI;
        adUpperBound[i] = (float) (i+1) / nI;
    }
    
    double adUniforms[nI];
    runif_vec(&nI, adLowerBound, adUpperBound, adUniforms);

    inverse_cdf_weights(nW, adWeights, nI, adUniforms, anIndices);
}


void systematic_resample(int nW, double *adWeights, int nI, int *anIndices)
{
    int i;
    double adUniforms[nI];
    GetRNGstate();
    adUniforms[0] = runif(0, (float) 1/ nI);
    PutRNGstate();
    for (i=1; i<nI; i++) adUniforms[i] =  adUniforms[i-1]+ (float) 1 / nI;

    inverse_cdf_weights(nW, adWeights, nI, adUniforms, anIndices);
}

void residual_resample(int nW, double *adWeights, int nI, int *anIndices,
                       int nResidualResampleFunction)
{

    // Particles are deterministically resampled floor(weights*nSamples) times
    int i, anDeterministicReps[nW], nDeterministicReps=0;
    double adExpectedSamples[nW];
    for (i=0; i<nW; i++) 
    {        
        // Expected samples 
        adExpectedSamples[i]    = adWeights[i]* nI;

        // Truncate to get deterministically resampled particles
        anDeterministicReps[i]  = adExpectedSamples[i];
     
        // Increment number of deterministic reps
        nDeterministicReps     += anDeterministicReps[i];

        // Remaining weight for use in random resampling
        adWeights[i]            = adExpectedSamples[i]-anDeterministicReps[i];
    }
    if (nDeterministicReps > nI) 
        REprintf("C: residual_resample: too many deterministic reps\n");
   
    rep2id(anDeterministicReps, nDeterministicReps, anIndices);


    // Particles are then randomly sampled with remaining weight
    nI -= nDeterministicReps;
    renormalize(nW, 0, adWeights);
    switch (nResidualResampleFunction) 
    {
        case 1:
            stratified_resample( nW, adWeights, nI, &anIndices[nDeterministicReps]);
            break;
        case 2:
            multinomial_resample(nW, adWeights, nI, &anIndices[nDeterministicReps]);
            break;
        case 3:
            systematic_resample( nW, adWeights, nI, &anIndices[nDeterministicReps]);
            break;
        default:
            REprintf("C: residual_resample: no match for residual resampling function\n");
    }
       
}



