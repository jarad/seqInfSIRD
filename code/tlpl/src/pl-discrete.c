/* 
* Functions for particle learning in discrete-time stochastic chemical kinetic models
*/

#include <assert.h>
#include <R.h>
#include <Rmath.h>
#include "utility.h"
#include "gillespie.h"
#include "pl-utility.h"
#include "pl-discrete.h"
#include "Sckm.h"
#include "SckmSwarm.h"
#include "SckmParticle.h"


/* Calculate the predictive likelihood */

void calc_log_pred_like_R(const int *anY, const double *dTau, 
                          int *nSpecies, int *nRxns, int *anPre, int *anPost, double *adlMult,
                          int *anX, double *probA, double *probB, double *rateA, double *rateB, 
                          double *prob, double *rate,
                          double *logPredLike)
{   
    Sckm *sckm = newSckm(*nSpecies, *nRxns, anPre, anPost, adlMult);        
    SckmParticle *part = newSckmParticle(sckm, anX, probA, probB, rateA, rateB, prob, rate);

    double adHazardPart[*nRxns];
    hazard_part(sckm, anX, adHazardPart);
    *logPredLike = calc_log_pred_like(anY, *dTau, sckm, part, adHazardPart);

    deleteSckm(sckm);
    deleteSckmParticle(part);
}

double calc_log_pred_like(const int *anY, double dTau, Sckm *sckm, SckmParticle *particle,
                          double *adHazardPart) 
{
    int nr = sckm->r;

    double adP2[nr], dLogPredLik=0;
    for (int i=0; i<nr; i++) {
        adP2[i] = 1/(1+particle->rateB[i]/(particle->prob[i]*adHazardPart[i]*dTau));
        if (adP2[i]>0)
            dLogPredLik += dnbinom(anY[i], particle->rateA[i], 1-adP2[i], 1);
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
    int i, nr=sckm->r;
    double adP[nr], adTheta[nr];
    GetRNGstate();
    for (i=0;i<nr;i++) 
    {
        adP[i] = rbeta(adHyper[i], adHyper[i+nr]);
        adTheta[i] = rgamma(adHyper[i+2*nr], adHyper[i+3*nr]);
    }
    PutRNGstate();

    double adHazardPart[nr], adHazard[nr];
    hazard(sckm, adTheta, anX, dTau, adHazardPart, adHazard);

    // Forward simulate system
    int anRxnCount[nr];
    *nSuccess = 1-cond_discrete_sim_step(sckm, adHazard, anY, adP, nWhileMax, anRxnCount, anX);

    suff_stat_update(nr, anRxnCount, anY, adHazardPart, adHyper);

    return 0;
}




void discrete_all_particle_update_R(int *nSpecies, int *nRxns, int *anPre, int *anPost, double *adlMult,
                                  const int *anY, const double *dTau,
                                  int *nParticles, int *nWhileMax,
                                  int *anX, double *adHyper, int *anSuccess) 
{
    Sckm *sckm = newSckm(*nSpecies, *nRxns, anPre, anPost, adlMult);        
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




void tlpl_R(
           /* Data */
           int *nObs,
           int *anY, 
           double *adTau, 

           /* sckm */
           int *nSpecies, 
           int *nRxns, 
           int *anPre, 
           int *anPost,
           double *adlMult,

           /* Particles */
           int *nParticles, 

           /* Auxiliary */
           int *nResamplingMethod,
           int *nNonuniformity,
           double *dThreshold,
           int *nVerbose,

           /* Outputs */
           int *anX,
           double *adProbA,
           double *adProbB,
           double *adRateA,
           double *adRateB
           )
{
    int i;
    Sckm *sckm = newSckm(*nSpecies, *nRxns, anPre, anPost, adlMult);

    SckmSwarm **swarm = newSckmSwarms(sckm, *nParticles, *nObs, 
                                       anX, adProbA, adProbB, adRateA, adRateB);
    
    tlpl(*nObs, anY, adTau,
         sckm, swarm, 
         *nResamplingMethod, *nNonuniformity, *dThreshold, *nVerbose);

    deleteSckmSwarms(swarm, *nObs);
    deleteSckm(sckm);
}

int tlpl(int nObs, int *anY, double *adTau,
         Sckm *sckm, SckmSwarm **swarm,
         int nResamplingMethod, int nNonuniformity, double dThreshold, int nVerbose)
{
    int nr = sckm->r, ns = sckm->s, np = swarm[0]-> nParticles; 

    double *hp      = (double *) malloc( np * nr * sizeof(double));
    double *prob    = (double *) malloc( np * nr * sizeof(double));
    double *rate    = (double *) malloc( np * nr * sizeof(double));
    double *weights = (double *) malloc( np *      sizeof(double));

    // Pointers 
    int *cY; cY = anY;
    double *cTau; cTau = adTau;
    
    double *cP, *cR, *cHP;
    SckmParticle *cPart;


    int i,j,k,l;
    for (i=0; i<nObs; i++) 
    {
        if (nVerbose) Rprintf("Time point %d, %3.0f%% completed.\n", i+1, (double) (i+1)/nObs*100);

        // Sample observation probability for all particles
        GetRNGstate();
        for (j=0; j< np; j++) 
        {
            cPart = swarm[i]->pParticle[j];
            cPart->prob = &prob[j * nr];
            for (k=0; k<nr; k++)
            {
                cPart->prob[k] = rbeta(cPart->probA[k], cPart->probB[k]);
            }
        }
        PutRNGstate();

        // Calculate hazard parts 
        for (j=0; j<np; j++) 
        {
            cPart = swarm[i]->pParticle[j];
            hazard_part(sckm, cPart->state, &hp[j*nr]);
            weights[j] = calc_log_pred_like(cY, *cTau, sckm, cPart, &hp[j*nr]);
        }
        

        // Update pointers
        cY += nr;
        cTau++; 
    }


    free(weights);
    free(rate);
    free(prob);
    free(hp);

    return 0;   
}





