#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <Rmath.h>

#include "Sckm.h"
#include "SckmParticle.h"
#include "SckmSwarm.h"


SckmSwarm *newSckmSwarm(Sckm *sckm, int _nParticles,
                        int *_state,
                        double *_probA, double *_probB, double *_rateA, double *_rateB,
                        double *_prob, double *_rate)
{
    SckmSwarm *swarm;
    swarm = (SckmSwarm *) malloc(sizeof(SckmSwarm));
    swarm->nParticles = _nParticles;
    swarm->nStates    = sckm->s;
    swarm->nRxns      = sckm->r;

    // Allocate dWeights and make uniform
    int i;
    double weight = 1.0/_nParticles;
    swarm->dWeights = (double *) malloc(_nParticles * sizeof(double));
    for (i=0; i< _nParticles; i++) swarm->dWeights[i] = weight;
    swarm->logWeights = 0;        // not log dWeights
    swarm->normalizedWeights = 1; // dWeights are normalized

    // Allocate and fill particles
    swarm->aParticles = (SckmParticle **) malloc(_nParticles * sizeof(SckmParticle *));
    for (i=0; i< _nParticles; i++) 
    {
        swarm->aParticles[i] = newSckmParticle(sckm, _state, _probA, _probB, _rateA, _rateB, _prob, _rate);
    }
}

void deleteSckmSwarm(SckmSwarm *swarm)
{
    free(swarm->dWeights);
    for (int i=0; i< swarm->nParticles; i++) deleteSckmParticle(swarm->aParticles[i]);
    free(swarm->aParticles);
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

