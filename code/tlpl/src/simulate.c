/* 
* Functions for discrete-time compartment models
*/

#include <R.h>
#include <Rmath.h>
#include "utility.h"
#include "simulate.h"

/* Calculates the part of the hazard other than the fixed parameter */
void hazard_part(const int *nSpecies, const int *nRxns, const int *anPre, // system specific arguments 
                 const int *anX,                                          // current system state
                 int *anHazardPart)                                       // return: (partial) system hazard
{
    int i, j;
    for (i=0; i<*nRxns; i++) 
    {
        anHazardPart[i] = 1;
        for (j=0; j<*nSpecies; j++) 
        {
            anHazardPart[i] *= choose(anX[j], anPre[i* *nSpecies+j]); 
        }    
    }       
}

/* Calculates the hazard for the next reaction */
void hazard(const int *nSpecies, const int *nRxns, const int *anPre, const double *adTheta,   
            const int *anX, 
            const double *dTau,  
            int *anHazardPart, double *adHazard)                   // return: hazard
{
    hazard_part(nSpecies, nRxns, anPre, anX, anHazardPart);
    int i;
    for (i=0; i<*nRxns; i++) 
    {
        adHazard[i] = adTheta[i]*anHazardPart[i]* *dTau;
    }
//    Rprintf("\n");
}


/* Simulates a set of reactions */
void rpois_vec(const int *nLength,                      
               const double *adMean,                     
               int *anPoisson)                         // return: number of reactions
{
    int i;
    GetRNGstate();
    for (i=0; i<*nLength; i++) anPoisson[i]=rpois(adMean[i]);
    PutRNGstate(); 
}


/* Updates the species according to the stoichiometry */
void update_species(const int *nSpecies, const int *nRxns,     
                    const int *anStoich, const int *anRxnCount,    
                    int *anX)                               // return: updated species
{
    int i,j;
    for (i=0; i<*nSpecies; i++) 
    {
        for (j=0; j<*nRxns; j++) 
        {
            anX[i] += anStoich[*nSpecies*j+i]*anRxnCount[j];
        }    
    } 
}

/* Forward simulate ahead one time-step */
void sim_one_step(const int *nSpecies, const int *nRxns, const int *anStoich, 
                  const double *adHazard,                 
                  const int *nWhileMax,                          
                  int *anRxnCount, int *anX)                                 // return: updated species
{
    int  whileCount=0, anTempX[*nSpecies];
    while (1) 
    {
        // Copy current state for temporary use
        copy(*nSpecies, anX, anTempX);

        // Get number of reactions
        rpois_vec(nRxns, adHazard, anRxnCount);

        // Update species
        update_species(nSpecies, nRxns, anStoich, anRxnCount, anTempX);

        // Test if update has any negative species
        if (!anyNegative(*nSpecies, anTempX)) 
        {
            // Copy temporary state into current state
            copy(*nSpecies, anTempX, anX);
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
    int i, nSO=0, anRxnCount[*nRxns], anHazardPart[*nRxns];
    double adHazard[*nRxns];
    for (i=0; i<*nSteps;i++)
    {
        // Calculate hazard based on current state
        hazard(nSpecies, nRxns, anPre, adTheta, &anX[nSO], &dTau[i], anHazardPart, adHazard);

        // Forward simulate the system
        sim_one_step(nSpecies, nRxns, anStoich, adHazard, nWhileMax, anRxnCount, &anX[nSO]);

        // Copy state for next step
        copy(*nSpecies, &anX[nSO], &anX[nSO+ *nSpecies]);
        nSO += *nSpecies;
    }
}


