#include "Sckm.h"

#ifndef SCKM_PARTICLE_H
#define SCKM_PARTICLE_H

typedef struct SckmParticle {
    int *state;
    double *probA, *probB, *rateA, *rateB; // hyperparameters 
    double *prob, *rate;                   // sampled parameters  
} SckmParticle;

SckmParticle *newSckmParticle(Sckm *sckm,
                              int *_state,
                              double *_probA, double *_probB, double *_rateA, double *_rateB,
                              double *_prob, double *_rate);

void setSckmParticle(SckmParticle *particle,
                     int *_state,
                     double *_probA, double *_probB, double *_rateA, double *_rateB,
                     double *_prob, double *_rate);

void deleteSckmParticle(SckmParticle *particle);


#endif // SCKM_PARTICLE_H
