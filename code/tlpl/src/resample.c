#include <R.h>
#include <Rmath.h>

/* Renormalizes the vector adWeights to sum to 1. If *nLog \ne 0, the adWeights are assumed to be log adWeights */
void renormalize(const int *n, const int *nLog, double *adWeights) 
{
  int i;

  if (*nLog) {
    double max=adWeights[0];
    for (i=1; i<*n; i++) if (adWeights[i]>max) max = adWeights[i]; 
    for (i=0; i<*n; i++) adWeights[i] = exp(adWeights[i]-max);
  }

  double sum=0;
  for (i=0; i<*n; i++) sum += adWeights[i];
  for (i=0; i<*n; i++) adWeights[i] /= sum;
}


/* Performs multinomial resampling on the adWeights. Returns the anIndices for nSamples samples */
void multinomial_resample(const int *n, const int *nSamples, const double *adWeights, int *anIndices) 
{
  int i, j;
  double cusum[*n];

  /* Calculate the cumulative sum of the adWeightseight series */
  cusum[0] = adWeights[0];
  for (i=1; i<*n; i++) cusum[i] = cusum[i-1]+adWeights[i];
  
  GetRNGstate();
  for (i=0; i<*nSamples; i++) {
    j=0;
    while (cusum[j]< runif(0,1)) { j++; }
    anIndices[i] = j;
  }
  PutRNGstate();
}

