/* 
* Functions for particle learning in discrete-time stochastic chemical kinetic models
*/

#include <R.h>
#include <Rmath.h>
#include "utility.h"
#include "gillespie.h"
#include "pl-utility.h"
#include "pl-discrete.h"
#include "Sckm.h"


/* Calculate the predictive likelihood */

void calc_log_pred_like_R(int *nSpecies, int *nRxns, int *anPre, int *anPost,
                             const int *anY, const int *anX, const double *adP, const double *adHyper,
                             double *logPredLike)
{   
    Sckm *sckm = newSckm(*nSpecies, *nRxns, anPre, anPost);        
    *logPredLike = calc_log_pred_like(sckm, anY, anX, adP, adHyper);
    deleteSckm(sckm);
}

double calc_log_pred_like(Sckm *sckm, const int *anY /* Y_{t+1} */, const int *anX /* X_t */,
                          const double *adP, const double *adHyper)   // only hyperparameters related to rates    
                                
{
    int nRxns = sckm->r;
    int anHazardPart[nRxns];
    hazard_part(sckm, anX, anHazardPart);

    double adP2[nRxns], dLogPredLik=0;
    for (int i=0; i<nRxns; i++) {
        adP2[i] = 1/(1+adHyper[i]/(adP[i]*anHazardPart[i]));
        dLogPredLik += dnbinom(adHyper[i+nRxns], anY[i], adP2[i], 1);
    }
    return dLogPredLik;
}


/* Sample from the state conditional on the observations */
int cond_discrete_sim_step(Sckm *sckm, const double *adHazard, const int *anY, const double *adP, 
                       int nWhileMax, int *anRxnCount, int *anX)
{
    // update hazard by probability of not observing
    int i, nRxns=sckm->r, nSpecies=sckm->s;
    double adHazardTemp[nRxns];
    for (i=0; i<nRxns; i++) adHazardTemp[i] = adHazard[i] * (1-adP[i]); 
    
    int whileCount=0, anTempX[nSpecies], anUnobservedRxnCount[nRxns], anTotalRxns[nRxns];
    while (1) 
    {
        memcpy(anTempX, anX, nSpecies*sizeof(int));

        // Sample unobserved reactions and add to observed reactions
        for (i=0; i<nRxns; i++) 
        {
            anUnobservedRxnCount[i] = rpois(adHazardTemp[i]);
            anTotalRxns[i] = anUnobservedRxnCount[i]+anY[i];
        }

        update_species(sckm, anTotalRxns, anTempX);

        if (!anyNegative(nSpecies, anTempX)) 
        {
            memcpy(anX, anTempX, nSpecies*sizeof(int));
            memcpy(anRxnCount, anTotalRxns, nRxns*sizeof(int));
            return 0;
        }

        // Limit how long the simulation tries to find a non-negative update
        whileCount++;
        if (whileCount>nWhileMax) 
            return 1;
            // error("C:cond_discrete_sim_step: Too many unsuccessful simulation iterations.");
    }
    return 0;
}  


/* Particle learning update for a single particle */
int discrete_particle_update(Sckm *sckm, const int *anY, double dTau, int nWhileMax,
                              int *anX, double *adHyper, int *nSuccess)
{
    // Sample parameters
    int i, nRxns=sckm->r;
    double adP[nRxns], adTheta[nRxns];
    GetRNGstate();
    for (i=0;i<nRxns;i++) 
    {
        adP[i] = rbeta(adHyper[i], adHyper[i+nRxns]);
        adTheta[i] = rgamma(adHyper[i+2*nRxns], adHyper[i+3*nRxns]);
    }
    PutRNGstate();

    int    anHazardPart[nRxns];
    double adHazard[    nRxns];
    hazard(sckm, adTheta, anX, dTau, anHazardPart, adHazard);

    // Forward simulate system
    int anRxnCount[nRxns];
    *nSuccess = 1-cond_discrete_sim_step(sckm, adHazard, anY, adP, nWhileMax, anRxnCount, anX);

    suff_stat_update(nRxns, anRxnCount, anY, anHazardPart, adHyper);

    return 0;
}




void discrete_all_particle_update_R(int *nSpecies, int *nRxns, int *anPre, int *anPost, 
                                  const int *anY, const double *dTau,
                                  int *nParticles, int *nWhileMax,
                                  int *anX, double *adHyper, int *anSuccess) 
{
    Sckm *sckm = newSckm(*nSpecies, *nRxns, anPre, anPost);        
    discrete_all_particle_update(sckm, anY,  *dTau, *nParticles, *nWhileMax,
                                 anX,  adHyper, anSuccess);
    deleteSckm(sckm);

}

/* Particle learning update for all particles */
int discrete_all_particle_update(Sckm *sckm, const int *anY, double dTau,
                                  int nParticles, int nWhileMax,
                                  int *anX, double *adHyper, int *anSuccess) 
{
    for (int i=0; i< nParticles; i++) 
    {
        discrete_particle_update(sckm, anY, dTau, nWhileMax, &anX[i* (sckm->s)], 
                                 &adHyper[i* 4*(sckm->r)], &anSuccess[i]); // 4 hyper parameters per reaction
    }
    return 0;
}





