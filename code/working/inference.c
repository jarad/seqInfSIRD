/* 
* Functions for discrete-time compartment models
*/

#include <R.h>
#include <Rmath.h>

/* Calculates the part of the hazard other than the fixed parameter */
void hazard_part(int *nSpecies, int *nRxns, int *anX, int *anPre, int *anHp) 
{
    int i, j;
    for (i=0; i<*nRxns; i++) 
    {
        anHp[i] = 1;
        for (j=0; j<*nSpecies; j++) 
        {
            anHp[i] *= choose(anX[j], anPre[i* *nSpecies+j]);            
        }    
    }       
}

/* Calculates the hazard for the next reaction */
void hazard(int *nSpecies, int *nRxns, int *anX, int *anPre, double *adTheta, int *anHp, double *adH) 
{
    hazard_part(nSpecies, nRxns, anX, anPre, anHp);
    int i;
    for (i=0; i<*nRxns; i++) adH[i] = adTheta[i]*anHp[i];
}


/* Simulates a set of reactions */
void sim_poisson(int *nRxns, double *adH, int *anRxns) 
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
        }    
    } 
}

/* Moves the system ahead one time-step */
void sim_one_step(int *nSpecies, int *nRxns, int *anX, int *anStoich, int *anRxns, double *adH) 
{
    int i, whileCount=0, anTempX[*nSpecies];
    while (1) 
    {
        // Copy current state for temporary use
        for (i=0; i<*nSpecies; i++) anTempX[i] = anX[i];

        // Get number of reactions
        sim_poisson(nRxns, adH, anRxns);

        // Update species
        update_species(nSpecies, nRxns, anTempX, anStoich, anRxns);

        // Test if update has any negative species
        if (!anyNegative(*nSpecies, anTempX)) 
        {
            // Copy temporary state into current state
            for (i=0; i<*nSpecies; i++) anX[i] = anTempX[i];
            break;
        }

        // Limit how long the simulation tries to find a non-negative update
        whileCount++;
        if (whileCount>1000) error("C:sim_one_step: Too many unsuccessful simulation iterations.");
    }
} 


void inf_one_step(int *nSpecies, int *nRxns, int *anX, int *anPre, int *anRxns, int *anY, int *anHp, int *anHyper) 
{
    
}



void one_step(int *nSpecies, int *nRxns, int *anX, int *anPre, int *anStoich, 
              double *adTheta, int *anY, int *anRxns, int *anHyper)
{
    // Calculate reaction hazard 
    int    anHp[*nRxns];
    double adH [*nRxns];
    hazard(nSpecies, nRxns, anX, anPre, adTheta, anHp, adH);

    // Forward simulate system
    sim_one_step(nSpecies, nRxns, anX, anStoich, anRxns, adH);

    // Inference for this simulated step
    inf_one_step(nSpecies, nRxns, anX, anPre, anRxns, anY, anHp, anHyper);
}


