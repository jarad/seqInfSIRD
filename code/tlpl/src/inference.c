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

/* Forward simulate ahead one time-step */
void sim_one_step(int *nSpecies, int *nRxns, int *anX, int *anStoich, int *anRxns, double *adH, 
                  int *nWhileMax) 
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
        if (whileCount>*nWhileMax) error("C:sim_one_step: Too many unsuccessful simulation iterations.");
    }
} 


/* Sample from the state conditional on the observations */
void cond_sim_one_step(int *nSpecies, int *nRxns, int *anX, int *anStoich, int *anRxns, 
                       double *adH, int *anY, double *adP, int *nWhileMax)
{
    int i, whileCount=0, anTempX[*nSpecies], anTotalRxns[*nRxns];
    for (i=0; i<*nRxns; i++) adH[i] *= (1-adP[i]); // update hazard by probability of not observing
    while (1) 
    {
        // Copy current state for temporary use
        for (i=0; i<*nSpecies; i++) anTempX[i] = anX[i];

        // Get number of reactions
        sim_poisson(nRxns, adH, anRxns);
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
void inf_one_step(int *nRxns, int *anRxns, int *anY, int *anHp, int *anHyper) 
{
    int i,j;
    int o1=*nRxns,o2=2**nRxns,o3=3**nRxns; // offsets
    for (i=0; i<*nRxns; i++) 
    {
        anHyper[i   ] += anY[i];           // Beta (alpha)        - observations
        anHyper[i+o1] += anRxns[i]-anY[i]; // Beta (beta)         - unobserved
        anHyper[i+o2] += anRxns[i];        // Gamma shape (alpha) - transitions
        anHyper[i+o3] += anHp[i];          // Gamma rate (beta)   - expected transitions
    }
}


/* Particle learning update for a single particle */
void one_step_single_particle(int *nSpecies, int *nRxns, int *anPre, int *anStoich, int *anY,
                              // Particle specific details
                              int *anX, double *adTheta, int *anRxns, int *anHyper, 
                              int *nWhileMax)
{
    // Calculate reaction hazard 
    int    anHp[*nRxns];
    double adH [*nRxns];
    hazard(nSpecies, nRxns, anX, anPre, adTheta, anHp, adH);

    // Forward simulate system
    sim_one_step(nSpecies, nRxns, anX, anStoich, anRxns, adH, nWhileMax);

    // Inference for this simulated step
    inf_one_step(nRxns, anRxns, anY, anHp, anHyper);
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

