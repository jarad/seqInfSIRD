#include <R.h>
#include <Rmath.h>
#include "pl-utility.h"

/* Particle learning sufficient statistic update */
int suff_stat_update(int nRxns, const int *anRxnCount, const int *anY, const double *anHazardPart, 
                      double *adHyper) 
{
    int i,j;
    int o1=nRxns,o2=2*nRxns,o3=3*nRxns; // offsets
    for (i=0; i<nRxns; i++) 
    {
        adHyper[i   ] += anY[i];               // Beta (alpha)        - observations
        adHyper[i+o1] += anRxnCount[i]-anY[i]; // Beta (beta)         - unobserved
        adHyper[i+o2] += anRxnCount[i];        // Gamma shape (alpha) - transitions
        adHyper[i+o3] += anHazardPart[i];      // Gamma rate (beta)   - expected transitions
    }
    return 0;
}



