/* 
* Functions for discrete-time compartment models
*/

#include <R.h>
#include <Rmath.h>

/* Calculates the part of the hazard other than the fixed parameter */
void hazard_part(int *nSpecies, int *nRxns, int *anX, int *anPre, int *anH) 
{
    int i, j, nOffset=0;
    for (i=0; i<*nRxns; i++) 
    {
        anH[i] = 1;
        for (j=0; j<*nSpecies; j++) 
        {
            anH[i] *= choose(anX[j], anPre[nOffset+j]);            
        }
        nOffset += *nSpecies;
    }       
}

/* Calculates the hazard for the next reaction */
void hazard(int *nSpecies, int *nRxns, int *anX, int *anPre, double *adTheta, double *adH) 
{
    int anHp[*nRxns];
    hazard_part(nSpecies, nRxns, anX, anPre, anHp);
    int i;
    for (i=0; i<*nRxns; i++) adH[i] = adTheta[i]*anHp[i];
}


/* Simulates a set of reactions */
void sim_poisson(int *nSpecies, int *nRxns, int *anX, int *anPre, double *adH, int *anRxns) 
{
    int i;
    GetRNGstate();
    for (i=0; i<*nRxns; i++) anRxns[i]=rpois(adH[i]);
    PutRNGstate(); 
}

int anyNegative(int n, int *v) 
{
    int i;
    for (i=0; i<n; i++) 
    {
        if (v[i]<0) return 1;
    }
    return 0;
}

void update_species(int *nSpecies, int *nRxns, int *anX, int *anStoich, int *anRxns) 
{
    int i,j;
    for (i=0; i<*nSpecies; i++) 
    {
        for (j=0; j<*nRxns; j++) 
        {
            anX[i] += anStoich[*nSpecies*j+i]*anRxns[j];
        }    } 
}

/* Moves the system ahead one time-step */
void sim_one_step(int *nSpecies, int *nRxns, int *anX, int *anPre, int *anStoich, 
              double *adTheta) 
{
    double adH[*nRxns];
    hazard(nSpecies, nRxns, anX, anPre, adTheta, adH);

    int i, whileCount=0, anRxns[*nRxns], anNewX[*nSpecies];
    while (1) 
    {
        for (i=0; i<*nSpecies; i++) anNewX[i] = anX[i];
        sim_poisson(nSpecies, nRxns, anX, anPre, adH, anRxns);
        update_species(nSpecies, nRxns, anNewX, anStoich, anRxns);
        if (!anyNegative(*nSpecies, anNewX)) 
        {
            for (i=0; i<*nSpecies; i++) anX[i] = anNewX[i];
            break;
        }
        whileCount++;
        if (whileCount>1000) error("C:sim_one_step: Too many unsuccessful simulation iterations.");
    }
} 

