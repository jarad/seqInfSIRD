/* 
* Functions for particle learning in discrete-time stochastic chemical kinetic models
*/

#include <R.h>
#include <Rmath.h>
#include "utility.h"
#include "simulate.h"

/* Sample from the state conditional on the observations */
void cond_sim_one_step(const int *nSpecies, const int *nRxns, const int *anStoich,  
                       const double *adHazard, const int *anY, const double *adP, const int *nWhileMax,
                       int *anRxns, int *anX)
{
    // update hazard by probability of not observing
    int i;
    double adHazardTemp[*nRxns];
    for (i=0; i<*nRxns; i++) adHazardTemp[i] = adHazard[i] * (1-adP[i]); 

    int whileCount=0, anTempX[*nSpecies], anTotalRxns[*nRxns];
    while (1) 
    {
        // Copy current state for temporary use
        copy(*nSpecies, anX, anTempX);

        // Get unobserved reactions and add to observed reactions
        rpois_vec(nRxns, adHazardTemp, anRxns);
        for (i=0; i<*nRxns; i++) anTotalRxns[i] = anRxns[i]+anY[i];

        // Temporarily update species according to temporary reactions
        update_species(nSpecies, nRxns, anTempX, anStoich, anTotalRxns);

        // Test if update has any negative species
        if (!anyNegative(*nSpecies, anTempX)) 
        {
            // Copy successful state and number of reactions back for returning from function
            copy(*nSpecies, anTempX,      anX);
            copy(*nRxns   , anTotalRxns, anRxns);
            break;
        }

        // Limit how long the simulation tries to find a non-negative update
        whileCount++;
        if (whileCount>*nWhileMax) error("C:sim_one_step: Too many unsuccessful simulation iterations.");
    }

}  

/* Particle learning sufficient statistic update */
void inf_one_step(const int *nRxns, const int *anRxns, const int *anY, const int *anHazardPart, 
                  int *anHyper) 
{
    int i,j;
    int o1=*nRxns,o2=2**nRxns,o3=3**nRxns; // offsets
    for (i=0; i<*nRxns; i++) 
    {
        anHyper[i   ] += anY[i];           // Beta (alpha)        - observations
        anHyper[i+o1] += anRxns[i]-anY[i]; // Beta (beta)         - unobserved
        anHyper[i+o2] += anRxns[i];        // Gamma shape (alpha) - transitions
        anHyper[i+o3] += anHazardPart[i];  // Gamma rate (beta)   - expected transitions
    }
}


/* Particle learning update for a single particle */
void one_step_single_particle(const int *nSpecies, const int *nRxns, const int *anPre, const int *anStoich, 
                              const int *anY, const double *dTau,
                              int *anX, double *adTheta, int *anRxns, int *anHyper, 
                              int *nWhileMax)
{
    // Calculate reaction hazard 
    int    anHazardPart[*nRxns];
    double adHazard[    *nRxns];
    hazard(nSpecies, nRxns, anPre, adTheta, anX, dTau, anHazardPart, adHazard);

    // Forward simulate system
    sim_one_step(nSpecies, nRxns, anStoich, adHazard, nWhileMax, anRxns, anX);

    // Inference for this simulated step
    inf_one_step(nRxns, anRxns, anY, anHazardPart, anHyper);
}


/* Particle learning update for all particles */
void one_step(const int *nSpecies, const int *nRxns, const int *anPre, const int *anStoich, 
              const int *anY, const double *dTau,
              const int *nParticles, const int *nWhileMax,
              int *anX, double *adTheta, int *anRxns, int *anHyper) 
              
{
    int i, j, nSO=0, nRO=0, anX0[*nSpecies];
    for (i=0; i< *nParticles; i++) 
    { 
        for (j=0; j< *nSpecies; j++) anX0[j] = anX[nSO+j]; 
        one_step_single_particle(nSpecies, nRxns, anPre, anStoich, anY, dTau,
                                 &anX[nSO], &adTheta[nRO], &anRxns[nRO], &anHyper[2*nRO],
                                 nWhileMax);
        nSO += *nSpecies; // species offset
        nRO += *nRxns;    // rxn offset
    }
}

