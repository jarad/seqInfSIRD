/* 
* Functions for particle learning in discrete-time stochastic chemical kinetic models
*/

#include <R.h>
#include <Rmath.h>
#include "utility.h"
#include "simulate.h"
#include "pl-utility.h"

/* Sample from the state conditional on the observations */
void cond_discrete_sim_step(const int *nSpecies, const int *nRxns, const int *anStoich,  
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
        if (whileCount>*nWhileMax) 
            error("C:sim_discrete_particle_update: Too many unsuccessful simulation iterations.");
    }

}  


/* Particle learning update for a single particle */
void discrete_particle_update(const int *nSpecies, const int *nRxns, const int *anPre, const int *anStoich, 
                              const int *anY, const double *dTau, const int *nWhileMax,
                              int *anX, double *adHyper)
{
    // Sample parameters
    double adP[*nRxns], adTheta[*nRxns];
    sample_p(    nRxns,  adHyper           , adP);
    sample_theta(nRxns, &adHyper[2* *nRxns], adTheta);

    // Calculate reaction hazard 
    int    anHazardPart[*nRxns];
    double adHazard[    *nRxns];
    hazard(nSpecies, nRxns, anPre, adTheta, anX, dTau, anHazardPart, adHazard);

    // Forward simulate system
    int anRxns[*nRxns];
    cond_discrete_sim_step(nSpecies, nRxns, anStoich, adHazard, anY, adP, nWhileMax, anRxns, anX);

    // Inference for this simulated step
    suff_stat_update(nRxns, anRxns, anY, anHazardPart, adHyper);
}


/* Particle learning update for all particles */
void discrete_all_particle_update(const int *nSpecies, const int *nRxns, const int *anPre, const int *anStoich, 
                                  const int *anY, const double *dTau,
                                  const int *nParticles, const int *nWhileMax,
                                  int *anX, double *adHyper) 
              
{
    int i, j, nSO=0, nRO=0, anX0[*nSpecies];
    for (i=0; i< *nParticles; i++) 
    { 
        for (j=0; j< *nSpecies; j++) anX0[j] = anX[nSO+j]; 
        discrete_particle_update(nSpecies, nRxns, anPre, anStoich, anY, dTau, nWhileMax,
                                 &anX[nSO], &adHyper[2*nRO]);
        nSO += *nSpecies; // species offset
        nRO += *nRxns;    // rxn offset
    }
}

