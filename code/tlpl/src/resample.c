#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include "utility.h"
#include "resample.h"

/*********************** Utility functions - to be moved *******************************/


// used in qsort and stolen from http://en.allexperts.com/q/C-1587/Qsort-function.htm
int compare_doubles (const void *X, const void *Y)
{
       double x = *((double *)X);
       double y = *((double *)Y);

       if (x > y)
       {
               return 1;
       }
       else
       {
               if (x < y)
               {
                       return -1;
               }
               else
               {
                       return 0;
               }
       }
}


void is_increasing_wrap(int *n, const double *v, int *returned) 
{
    *returned = is_increasing(*n, v); 
}

int is_increasing(int n, const double *v) 
{
    int i; 
    for (i=1; i<n; i++)
    {
        if (v[i]<v[i-1]) return 0;
    }
    return 1;
}


void cumulative_sum_wrap(int *n, double *v) 
{
    cumulative_sum(*n, v);
}

int cumulative_sum(int n, double *v) 
{
    int i;
    for (i=1; i<n; i++) v[i] += v[i-1];
    return 0;
}


void rep2id_wrap(int *rep, int *sum, int *id) 
{
    rep2id(rep, *sum, id);
}

int rep2id(int *rep, int sum, int *id)
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

    return 0;
}


void inverse_cdf_weights_wrap(int *nW, double *adWeights, int *nU, double *adUniforms, int *anIndices)
{
    inverse_cdf_weights(*nW, adWeights, *nU, adUniforms, anIndices);
}


int inverse_cdf_weights(int nW, 
                         double *adWeights, 
                         int nU, 
                         double *adUniforms,
                         int *anIndices)
{
    if (!is_increasing(nU, adUniforms))    
        qsort(adUniforms, nU, sizeof(double), compare_doubles);

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

    return 0;   
}


/************************ Effective sample size functions ******************************/


void effective_sample_size_wrap(int *n, double *weights, double *returned)
{
    *returned = effective_sample_size(*n, weights);
}

double effective_sample_size(int n, double *weights)
{
    int i;
    double sum=0;
    for (i=0; i<n; i++) sum += weights[i]*weights[i];
    return 1/sum;
}



void coefficient_of_variation_wrap(int *n, double *weights, double *returned)
{
    *returned = coefficient_of_variation(*n, weights);
}

double coefficient_of_variation(int n, double *weights) 
{
    int i;
    double mean=0, var=0;
    for (i=0; i<n; i++) mean += weights[i];
    mean /= n;
    for (i=0; i<n; i++) var += R_pow_di(weights[i]-mean,2);
    return var/R_pow_di(mean,2); // why this as opposed to the sqrt of this?
}




void entropy_wrap(int *n, double *weights, double *returned) {
    *returned = entropy(*n, weights);
}

double entropy(int n, double *weights)
{
    int i;
    double sum;
    for (i=0; i<n; i++) sum += weights[i]*log2(weights[i]); // should add smallest constant within log2()
    return -sum;
}

/************************ Resampling functions *****************************************/




/* Renormalizes the vector w to sum to 1. 
 * If log \ne 0, the weights are assumed to be log weights */

void renormalize_wrap(int *n, int *log, double *w)
{
    renormalize(*n, *log, w);
}

int renormalize(int n, int log, double *w) 
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

    return 0;
}


void multinomial_resample_wrap( int *nW, double *adWeights, int *nI, int *anIndices)
{
    multinomial_resample(*nW, adWeights, *nI, anIndices);
}

int multinomial_resample(int nW, double *adWeights, int nI, int *anIndices) 
{
    REprintf("C: multinomial_resample: currently not working\n");
    int i, j;
    double cusum[nW], adUniforms[nI];

    GetRNGstate();
    for (i=0; i<nI; i++) adUniforms[i] = runif(0,1);
    PutRNGstate();

    inverse_cdf_weights(nW, adWeights, nI, adUniforms, anIndices);

    return 0;
}




void stratified_resample_wrap(int *nW, double *adWeights, int *nI, int *anIndices)
{
    stratified_resample(*nW, adWeights, *nI, anIndices);
}



int stratified_resample(int nW, double *adWeights, int nI, int *anIndices)
{
    int i;
    double adLowerBound[nI], adUpperBound[nI];
    for (i=0;i<nI;i++)
    {
        adLowerBound[i] = (float)  i    / nI;
        adUpperBound[i] = (float) (i+1) / nI;
    }
    
    double adUniforms[nI];
    runif_vec(nI, adLowerBound, adUpperBound, adUniforms);

    inverse_cdf_weights(nW, adWeights, nI, adUniforms, anIndices);

    return 0;
}




void systematic_resample_wrap(int *nW, double *adWeights, int *nI, int *anIndices)
{
    systematic_resample(*nW, adWeights, *nI, anIndices);
}


int systematic_resample(int nW, double *adWeights, int nI, int *anIndices)
{
    int i;
    double adUniforms[nI];
    GetRNGstate();
    adUniforms[0] = runif(0, (float) 1/ nI);
    PutRNGstate();
    for (i=1; i<nI; i++) adUniforms[i] =  adUniforms[i-1]+ (float) 1 / nI;

    inverse_cdf_weights(nW, adWeights, nI, adUniforms, anIndices);

    return 0;
}





void residual_resample_wrap(int *nW, double *adWeights, int *nI, int *anIndices, 
                            int *nResidualResampleFunction)
{
    residual_resample(*nW, adWeights, *nI, anIndices, *nResidualResampleFunction);
}

int residual_resample(int nW, double *adWeights, int nI, int *anIndices,
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
       
    return 0;
}




