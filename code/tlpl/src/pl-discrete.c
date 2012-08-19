/* 
* Functions for particle learning in discrete-time stochastic chemical kinetic models
*/

#include <R.h>
#include <Rmath.h>
#include "utility.h"
#include "simulate.h"
#include "pl-utility.h"
#include "pl-discrete.h"


/* Calculate the predictive likelihood */
int calculate_log_predictive_likelihood(int nSpecies, int nRxns, const int *anPre,
                                    const int *anY, // Y_{t+1}
                                    const int *anX, // X_t
                                    const double *adP, 
                                    const double *adHyper,       // Hyperparameters for rates 
                                    double *adLogPredLike)       // Return
                                
{
    int i, anHazardPart[nRxns];
    hazard_part(nSpecies, nRxns, anPre, anX, anHazardPart);

    double adP2[nRxns];
    for (i=0; i<nRxns; i++) {
        adP2[i] = 1/(1+adHyper[i]/(adP[i]*anHazardPart[i]));
        adLogPredLike[i] = dnbinom(adHyper[i+nRxns], anY[i], adP2[i], 1);
    }
}


/* Sample from the state conditional on the observations */
int cond_discrete_sim_step(int nSpecies, int nRxns, const int *anStoich,  
                       const double *adHazard, const int *anY, const double *adP, int nWhileMax,
                       int *anRxnCount, int *anX)
{
    // update hazard by probability of not observing
    int i;
    double adHazardTemp[nRxns];
    for (i=0; i<nRxns; i++) adHazardTemp[i] = adHazard[i] * (1-adP[i]); 
    
    int whileCount=0, anTempX[nSpecies], anUnobservedRxnCount[nRxns], anTotalRxns[nRxns];
    while (1) 
    {
        // Copy current state for temporary use
        copy(nSpecies, anX, anTempX);

        // Get unobserved reactions and add to observed reactions
        rpois_vec(nRxns, adHazardTemp, anUnobservedRxnCount);
        for (i=0; i<nRxns; i++) anTotalRxns[i] = anUnobservedRxnCount[i]+anY[i];

        // Temporarily update species according to temporary reactions
        update_species(nSpecies, nRxns, anStoich, anTotalRxns, anTempX);

        // Test if update has any negative species
        if (!anyNegative(nSpecies, anTempX)) 
        {
            // Copy successful state and number of reactions back for returning from function
            copy(nSpecies, anTempX,     anX);
            copy(nRxns   , anTotalRxns, anRxnCount);
            return 1;
        }

        // Limit how long the simulation tries to find a non-negative update
        whileCount++;
        if (whileCount>nWhileMax) 
            return 0;
            // error("C:cond_discrete_sim_step: Too many unsuccessful simulation iterations.");
    }

}  


/* Particle learning update for a single particle */
void discrete_particle_update(int nSpecies, int nRxns, const int *anPre, const int *anStoich, 
                              const int *anY, double dTau, int nWhileMax,
                              int *anX, double *adHyper, int *nSuccess)
{
    // Sample parameters
    double adP[nRxns], adTheta[nRxns];
    rbeta_vec( nRxns,  adHyper,          &adHyper[  nRxns], adP);     // sample probabilities
    rgamma_vec(nRxns, &adHyper[2*nRxns], &adHyper[3*nRxns], adTheta); // sample reaction rates

    int i;
    Rprintf("  Sampled parameters: \n");
    Rprintf("    P: ");
    for (i=0; i<nRxns; i++) Rprintf("%2.8f ", adP[i]);
    Rprintf("\n    Theta: ");
    for (i=0; i<nRxns; i++) Rprintf("%2.8f ", adTheta[i]);
    Rprintf("\n");    

    // Calculate reaction hazard 
    int    anHazardPart[nRxns];
    double adHazard[    nRxns];
    hazard(nSpecies, nRxns, anPre, adTheta, anX, dTau, anHazardPart, adHazard);

    // Forward simulate system
    int anRxnCount[nRxns];
    *nSuccess = cond_discrete_sim_step(nSpecies, nRxns, anStoich, adHazard, anY, adP, nWhileMax, anRxnCount, anX);

    // Inference for this simulated step
    suff_stat_update(nRxns, anRxnCount, anY, anHazardPart, adHyper);
}






void discrete_all_particle_update_wrap(int *nSpecies, int *nRxns, const int *anPre, const int *anStoich, 
                                  const int *anY, const double *dTau,
                                  int *nParticles, int *nWhileMax,
                                  int *anX, double *adHyper, int *anSuccess) 
{
    discrete_all_particle_update(*nSpecies, *nRxns, anPre, anStoich, 
                                  anY,  *dTau,
                                  *nParticles, *nWhileMax,
                                  anX,  adHyper, anSuccess);

}

/* Particle learning update for all particles */
void discrete_all_particle_update(int nSpecies, int nRxns, const int *anPre, const int *anStoich, 
                                  const int *anY, double dTau,
                                  int nParticles, int nWhileMax,
                                  int *anX, double *adHyper, int *anSuccess) 
{
    int i;
    for (i=0; i< nParticles; i++) 
    {
        Rprintf("Particle %d:\n", i);
        discrete_particle_update(nSpecies, nRxns, anPre, anStoich, anY, dTau, nWhileMax,
                                 &anX[i* nSpecies], 
                                 &adHyper[i* 4*nRxns], &anSuccess[i]); // 4 hyper parameters per reaction
    }
}





