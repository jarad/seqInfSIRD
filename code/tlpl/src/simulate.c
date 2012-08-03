/* 
* Functions for discrete-time compartment models
*/

#include <R.h>
#include <Rmath.h>
#include "utility.h"

/* Calculates the part of the hazard other than the fixed parameter */
void hazard_part(int *nSpecies, int *nRxns, int *anPre, // these should be const
                 int *anX,                              // and this
                 int *anHp)                             // return hazard part
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
void hazard(int *nSpecies, int *nRxns, int *anPre,    // these should be const
            int *anX, double *adTheta, double *dTau,  // and these
            int *anHp, double *adH)                   // return hazard
{
    hazard_part(nSpecies, nRxns, anPre, anX, anHp);
    int i;
    for (i=0; i<*nRxns; i++) adH[i] = adTheta[i]*anHp[i]* *dTau;
}


/* Simulates a set of reactions */
void sim_poisson(int *nRxns,                       // this should be const
                 double *adH,                      // and this
                 int *anRxns)                      // return number of reactions
{
    int i;
    GetRNGstate();
    for (i=0; i<*nRxns; i++) anRxns[i]=rpois(adH[i]);
    PutRNGstate(); 
}


/* Updates the species according to the stoichiometry */
void update_species(int *nSpecies, int *nRxns,     // these should be const
                    int *anStoich, int *anRxns,    // and these
                    int *anX)                      // return updated species
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

/* Forward simulate ahead one time-step */
void sim_one_step(int *nSpecies, int *nRxns, int *anStoich, // const
                  int *anRxns, double *adH,                 // const
                  int *nWhileMax,                           // const
                  int *anX)                                 // return updated species
{
    int  whileCount=0, anTempX[*nSpecies];
    while (1) 
    {
        // Copy current state for temporary use
        copy_int(*nSpecies, anX, anTempX);

        // Get number of reactions
        sim_poisson(nRxns, adH, anRxns);

        // Update species
        update_species(nSpecies, nRxns, anStoich, anRxns, anTempX);

        // Test if update has any negative species
        if (!anyNegative(*nSpecies, anTempX)) 
        {
            // Copy temporary state into current state
            copy_int(*nSpecies, anTempX, anX);
            break;
        }

        // Limit how long the simulation tries to find a non-negative update
        whileCount++;
        if (whileCount>*nWhileMax) error("C:sim_one_step: Too many unsuccessful simulation iterations.");
    }
} 


void sim(int *nSpecies, int *nRxns, int *anStoich, int *anPre, double *adTheta,
         double *dTau, int *nSteps, 
         int *nWhileMax,
         int *anX)
{
    int i, nSO=0, anRxns[*nRxns], anHp[*nRxns];
    double adH[*nRxns];
    for (i=0; i<*nSteps;i++)
    { 
        // Calculate hazard based on current state
        hazard(nSpecies, nRxns, anPre, &anX[nSO], adTheta, dTau, anHp, adH);

        // Forward simulate the system
        sim_one_step(nSpecies, nRxns, anStoich, anRxns, adH, nWhileMax, &anX[nSO]);

        // Copy state for next step
        copy_int(*nSpecies, &anX[nSO], &anX[nSO+ *nSpecies]);
        nSO += *nSpecies;
    }
}


