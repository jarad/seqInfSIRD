#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Sckm.h"
#include "SckmParticle.h"

SckmParticle *newSckmParticle(Sckm *sckm,
                              int *_state,
                              double *_probA, double *_probB, double *_rateA, double *_rateB,
                              double *_prob, double *_rate) 
{
    int nSpecies = sckm->s, nRxns = sckm->r;
    size_t intv = nSpecies*sizeof(int), douv = nRxns*sizeof(double);

    SckmParticle *particle;
    particle        = (SckmParticle *) malloc(sizeof(SckmParticle));
    particle->state = (int *)          malloc(intv);
    particle->probA = (double *)       malloc(douv);
    particle->probB = (double *)       malloc(douv);
    particle->rateA = (double *)       malloc(douv);
    particle->rateB = (double *)       malloc(douv);
    particle->prob  = (double *)       malloc(douv);
    particle->rate  = (double *)       malloc(douv);

    memcpy(particle->state, _state, intv);
    memcpy(particle->probA, _probA, douv);
    memcpy(particle->probB, _probB, douv);
    memcpy(particle->rateA, _rateA, douv);
    memcpy(particle->rateB, _rateB, douv);
    memcpy(particle->prob , _prob , douv);
    memcpy(particle->rate , _rate , douv);

    return(particle);
}


// Switches pointers to point to a new particle
void setSckmParticle(SckmParticle *particle,
                     int *_state,
                     double *_probA, double *_probB, double *_rateA, double *_rateB,
                     double *_prob, double *_rate) 
{
    particle->state = _state;
    particle->probA = _probA;
    particle->probB = _probB;
    particle->rateA = _rateA;
    particle->rateB = _rateB;
    particle->prob  = _prob;
    particle->rate  = _rate;
}





void deleteSckmParticle(SckmParticle *particle)
{
    assert(particle);
    assert(particle->state); free(particle->state);
    assert(particle->probA); free(particle->probA);
    assert(particle->probB); free(particle->probB);
    assert(particle->rateA); free(particle->rateA);
    assert(particle->rateB); free(particle->rateB);
    assert(particle->prob ); free(particle->prob );
    assert(particle->rate ); free(particle->rate );
    free(particle);
}

