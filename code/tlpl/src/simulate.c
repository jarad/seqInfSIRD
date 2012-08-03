/* 
* Functions for discrete-time compartment models
*/

#include <R.h>
#include <Rmath.h>
#include "utility.h"

/* Calculates the part of the hazard other than the fixed parameter */
void hazard_part(const int *nSpecies, const int *nRxns, const int *anPre, // system specific arguments 
                 const int *anX,                                          // current system state
                 int *anHp)                                               // return: (partial) system hazard
{
    int i, j;
    for (i=0; i<*nRxns; i++) 
    {
        anHp[i] = 1;
        for (j=0; j<*nSpecies; j++) 
        {
            anHp[i] *= choose(anX[j], anPre[i* *nSpecies+j]); 
            Rprintf("%d %d: %d %d %d\n", i, j, anX[j], anPre[i* *nSpecies+j], anHp[i]);           
        }    
    }       
}

/* Calculates the hazard for the next reaction */
void hazard(const int *nSpecies, const int *nRxns, const int *anPre, const double *adTheta,   
            const int *anX, 
            const double *dTau,  
            int *anHp, double *adH)                   // return: hazard
{
    hazard_part(nSpecies, nRxns, anPre, anX, anHp);
    int i;
    for (i=0; i<*nRxns; i++) 
    {
        adH[i] = adTheta[i]*anHp[i]* *dTau;
        Rprintf("%d: %f %d %f %f\n", i, adTheta[i], anHp[i], *dTau, adH[i]);
    }
    Rprintf("\n");
}


/* Simulates a set of reactions */
void sim_poisson(const int *nRxns,                      
                 const double *adH,                     
                 int *anRxns)                         // return: number of reactions
{
    int i;
    GetRNGstate();
    for (i=0; i<*nRxns; i++) anRxns[i]=rpois(adH[i]);
    PutRNGstate(); 
}


/* Updates the species according to the stoichiometry */
void update_species(const int *nSpecies, const int *nRxns,     
                    const int *anStoich, const int *anRxns,    
                    int *anX)                               // return: updated species
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
void sim_one_step(const int *nSpecies, const int *nRxns, const int *anStoich, 
                  const double *adH,                 
                  const int *nWhileMax,                          
                  int *anRxns, int *anX)                                 // return: updated species
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


void sim(const int *nSpecies, const int *nRxns, const int *anStoich, const int *anPre, const double *adTheta,
         const double *dTau, const int *nSteps, 
         const int *nWhileMax,
         int *anX)
{
    int i, nSO=0, anRxns[*nRxns], anHp[*nRxns];
    double adH[*nRxns];
    for (i=0; i<*nSteps;i++)
    {
        Rprintf("Step %d:\n", i); 
        // Calculate hazard based on current state
        hazard(nSpecies, nRxns, anPre, adTheta, &anX[nSO], dTau, anHp, adH);

        // Forward simulate the system
        sim_one_step(nSpecies, nRxns, anStoich, adH, nWhileMax, anRxns, &anX[nSO]);

        // Copy state for next step
        copy_int(*nSpecies, &anX[nSO], &anX[nSO+ *nSpecies]);
        nSO += *nSpecies;
        Rprintf("\n\n\n");
    }
}


