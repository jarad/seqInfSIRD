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
#include "resample.h"


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
    int anUnobservedTransitions[nr], anTotalTransitions[nr];

    double rate[nr], mn[nr];
    double *hp   = (double *) malloc( np * nr * sizeof(double));
    double *prob = (double *) malloc( np * nr * sizeof(double));

    // Pointers 
    int *cY;      cY   = anY;   // current observations
    double *cTau; cTau = adTau; // current tau
    
    double *w;           // swarm weights
    SckmSwarm *s;        // current swarm
    SckmParticle *cPart; // current particle
    SckmParticle *nPart; // next particles


    int i,j,k,l, nDoResample, anResampledIndices[np], nAnyNegative;
    for (i=0; i<nObs; i++) 
    {
        if (nVerbose) Rprintf("Time point %d, %3.0f%% completed.\n", i+1, (double) (i+1)/nObs*100);

        // Update pointers
        s = swarm[i];
        w = s->adWeights; 

        // Sample observation probability for all particles
        GetRNGstate();
        for (j=0; j< np; j++) 
        {
            cPart = s->pParticle[j];
            cPart->prob = &prob[j * nr];
            for (k=0; k<nr; k++)
            {
                cPart->prob[k] = rbeta(cPart->probA[k], cPart->probB[k]);
            }
        }
        PutRNGstate();

        // Calculate particle weights         
        for (j=0; j<np; j++) 
        {
            cPart = s->pParticle[j];
            hazard_part(sckm, cPart->state, &hp[j*nr]);
            w[j] = calc_log_pred_like(cY, *cTau, sckm, cPart, &hp[j*nr]);
        }
        s->logWeights = 1;
        s->normalizedWeights = 0;
        
        // Resampling
        renormalize(s);
        nDoResample = doResample(np, w, nNonuniformity, dThreshold);
        if (nDoResample) {
            if (nVerbose>1) Rprintf(" Resampling.");
            resample(np, w, np, anResampledIndices, nResamplingMethod);
        } else {
            for (j=0; j<np; j++) anResampledIndices[j] = j;
        }

        // Update particles
        for (j=0; j<np; j++) 
        {
            if (nVerbose>1) Rprintf(" Particle %d\n", j);
            nPart = swarm[i+1]->pParticle[j];

            nAnyNegative = 1;
            while(nAnyNegative)
            {
                // The first time through take the resampled particle
                // on any other pass take a new particle
                k = nAnyNegative==1 ? anResampledIndices[j] : one_multinomial_sample(np, w);
                cPart = s->pParticle[k];               

                GetRNGstate();
                for (l=0; l<nr; l++)
                {
                    // Calculate expected transitions
                    rate[l] = rgamma(cPart->rateA[l], 1.0/cPart->rateB[l]); 
                    mn[l]   = (1-cPart->prob[l])*rate[l]*hp[k*nr+l];
                }

                // Two loops to match R
                for (l=0; l<nr; l++)
                {
                    // Sample transitions and update state
                    anUnobservedTransitions[l] = rpois(mn[l]);
                    anTotalTransitions[l] = anUnobservedTransitions[l]+cY[l];
                }
                PutRNGstate();

                memcpy(nPart->state, cPart->state, ns*sizeof(double));
                update_species(sckm, anTotalTransitions, nPart->state);

                // Check for negative states 
                for (l=0; l<ns; l++) 
                {
                    if (nPart->state[l] < 0) 
                    {
                        nAnyNegative++; 
                        nAnyNegative =0; // temporary
                        break;
                    }
                    nAnyNegative = 0;
                }
            } // while(nAnyNegative)

            // Update sufficient statistics
            for (l=0; l<nr; l++) 
            {
                nPart->probA[l] = cPart->probA[l] + cY[l];
                nPart->probB[l] = cPart->probB[l] + anUnobservedTransitions[l];
                nPart->rateA[l] = cPart->rateA[l] + anTotalTransitions[l];
                nPart->rateB[l] = cPart->rateB[l] + hp[k*nr+l];
            }
        } // Updated particles

        // Update pointers
        cY += nr;
        cTau++; 
    }

    free(prob);
    free(hp);

    return 0;   
}





