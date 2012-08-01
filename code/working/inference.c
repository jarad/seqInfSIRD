/* 
* Functions for discrete-time compartment models
*/

#include <R.h>

void hazard_part(int *nSpecies, int *nRxns, int *anX, int *anPre, int *anH) 
{
    int i, j, rOffset=0;
    for (i=0; i<*nRxns; i++) 
    {
        anH[i] = 1;
        for (j=0; j<*nSpecies; j++) 
        {
            switch (anPre[rOffset+j])
            {
                case 0:
                    break;
                case 1:
                    anH[i] *= anX[j];
                    break;
                case 2:
                    anH[i] *= anX[j]*(anX[j]-1)/2;
                    break;
                default:
                    error("Pre must have elements 0, 1, and 2.");
                    break;
             }
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




