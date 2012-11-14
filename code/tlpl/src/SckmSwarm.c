#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <R.h>
#include <Rmath.h>

#include "Sckm.h"
#include "SckmParticle.h"
#include "SckmSwarm.h"


SckmSwarm *newSckmSwarm(Sckm *sckm, int _nParticles,
                        int *_state,
                        double *_probA, double *_probB, double *_rateA, double *_rateB)
{
    SckmSwarm *swarm;
    swarm = (SckmSwarm *) malloc(sizeof(SckmSwarm));
    swarm->nParticles = _nParticles;
    swarm->nStates    = sckm->s;
    swarm->nRxns      = sckm->r;

    // Allocate dWeights and make uniform
    swarm->dWeights = (double *) malloc(_nParticles * sizeof(double));
    memset(swarm->dWeights, 0, _nParticles);

    swarm->logWeights = 1;        // log dWeights
    swarm->normalizedWeights = 1; // dWeights are normalized

    // Associate particle pointers
    swarm->pParticle = (SckmParticle *) malloc(sizeof(SckmParticle *));
    swarm->pParticle->state = _state;
    swarm->pParticle->probA = _probA;
    swarm->pParticle->probB = _probB;
    swarm->pParticle->rateA = _rateA;
    swarm->pParticle->rateB = _rateB;
}

void deleteSckmSwarm(SckmSwarm *swarm)
{
    // These free's created an error. But without them don't I have a memory leak?
    //free(swarm->pParticle);
    //free(swarm->dWeights);
    free(swarm);
}

int renormalize(SckmSwarm *swarm)
{
    if (swarm->normalizedWeights) return 0;

    int i, n=swarm->nParticles;
    double *w;
    w = swarm->dWeights; 

    if (swarm->logWeights) 
    {
        double max=w[0];
        for (i=1; i<n; i++) max = fmax2(max, w[i]);
        for (i=0; i<n; i++) w[i] = exp(w[i]-max);
    }
    swarm -> logWeights = 0;

    double sum=0;
    for (i=0; i<n; i++) sum += w[i];
    for (i=0; i<n; i++) w[i] /= sum;

    return 0;
}

