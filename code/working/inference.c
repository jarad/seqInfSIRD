/* 
* Functions for discrete-time compartment models
*/

#include <R.h>
#include <Rmath.h>

void hazard_part(int *nSpecies, int *nRxns, int *anX, int *anPre, int *anH) 
{
    int i, j, rOffset=0;
    for (i=0; i<*nRxns; i++) 
    {
        anH[i] = 1;
        for (j=0; j<*nSpecies; j++) 
        {
            anH[i] *= choose(anX[j], anPre[rOffset+j]);            
        }
        rOffset += *nSpecies;
    }       
}

void hazard(int *nSpecies, int *nRxns, int *anX, int *anPre, double *adH, double *adTheta) 
{
    int anHp[*nRxns];
    hazard_part(nSpecies, nRxns, anX, anPre, anHp);
    int i;
    for (i=0; i<*nRxns; i++) adH[i] = adTheta[i]*anHp[i];
}




