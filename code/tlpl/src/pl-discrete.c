/* 
* Functions for particle learning in discrete-time stochastic chemical kinetic models
*/

#include <R.h>
#include <Rmath.h>


/* Sample from the state conditional on the observations */
void cond_sim_one_step(int *nSpecies, int *nRxns, int *anStoich, int *anRxns, 
                       double *adHazard, int *anY, double *adP, int *nWhileMax,
                       int *anX)
{
    int i, whileCount=0, anTempX[*nSpecies], anTotalRxns[*nRxns];
    for (i=0; i<*nRxns; i++) adHazard[i] *= (1-adP[i]); // update hazard by probability of not observing
    while (1) 
    {
        // Copy current state for temporary use
        for (i=0; i<*nSpecies; i++) anTempX[i] = anX[i];

        // Get number of reactions
        sim_poisson(nRxns, adHazard, anRxns);
        for (i=0; i<*nRxns; i++) anTotalRxns[i] = anRxns[i]+anY[i];

        // Update species
        update_species(nSpecies, nRxns, anTempX, anStoich, anTotalRxns);

        // Test if update has any negative species
        if (!anyNegative(*nSpecies, anTempX)) 
        {
            // Copy temporary state into current state
            for (i=0; i<*nSpecies; i++) anX[i]    = anTempX[i];
            for (i=0; i<*nRxns   ; i++) anRxns[i] = anTotalRxns[i];
            break;
        }

        // Limit how long the simulation tries to find a non-negative update
        whileCount++;
        if (whileCount>*nWhileMax) error("C:sim_one_step: Too many unsuccessful simulation iterations.");
    }

}  

/* Particle learning sufficient statistic update */
void inf_one_step(int *nRxns, int *anRxns, int *anY, int *anHazardPart, int *anHyper) 
{
    int i,j;
    int o1=*nRxns,o2=2**nRxns,o3=3**nRxns; // offsets
    for (i=0; i<*nRxns; i++) 
    {
        anHyper[i   ] += anY[i];           // Beta (alpha)        - observations
        anHyper[i+o1] += anRxns[i]-anY[i]; // Beta (beta)         - unobserved
        anHyper[i+o2] += anRxns[i];        // Gamma shape (alpha) - transitions
        anHyper[i+o3] += anHazardPart[i];          // Gamma rate (beta)   - expected transitions
    }
}


/* Particle learning update for a single particle */
void one_step_single_particle(int *nSpecies, int *nRxns, int *anPre, int *anStoich, int *anY,
                              // Particle specific details
                              int *anX, double *adTheta, int *anRxns, int *anHyper, 
                              int *nWhileMax)
{
    // Calculate reaction hazard 
    int    anHazardPart[*nRxns];
    double adHazard [*nRxns];
    hazard(nSpecies, nRxns, anX, anPre, adTheta, anHazardPart, adHazard);

    // Forward simulate system
    sim_one_step(nSpecies, nRxns, anX, anStoich, anRxns, adHazard, nWhileMax);

    // Inference for this simulated step
    inf_one_step(nRxns, anRxns, anY, anHazardPart, anHyper);
}


/* Particle learning update for all particles */
void one_step(int *nSpecies, int *nRxns, int *anPre, int *anStoich, int *anY,
              // Particle specific arguments
              int *anX, double *adTheta, int *anRxns, int *anHyper, 
              int *nParticles, int *nWhileMax)
{
    int i, j, nSO=0, nRO=0, anX0[*nSpecies];
    for (i=0; i< *nParticles; i++) 
    { 
        for (j=0; j< *nSpecies; j++) anX0[j] = anX[nSO+j]; 
        one_step_single_particle(nSpecies, nRxns, anPre, anStoich, anY,
                                 &anX[nSO], &adTheta[nRO], &anRxns[nRO], &anHyper[2*nRO],
                                 nWhileMax);
        nSO += *nSpecies; // species offset
        nRO += *nRxns;     // rxn offset
    }
}

