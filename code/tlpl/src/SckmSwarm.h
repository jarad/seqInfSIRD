#ifndef SCKM_SWARM_H
#define SCKM_SWARM_H

#include "Sckm.h"
#include "SckmParticle.h"

typedef struct StochasticChemicalKineticModelSwarm {
    int nParticles, nStates, nRxns, logWeights, normalizedWeights;
    double *adWeights;
    SckmParticle **pParticle;
} SckmSwarm;


SckmSwarm **newSckmSwarms(Sckm *sckm, int _nParticles, int _nObs,
                        int *_state,
                        double *_probA, double *_probB, double *_rateA, double *_rateB);
void **deleteSckmSwarms(SckmSwarm **swarms, int nObs);

SckmSwarm *newSckmSwarm(Sckm *sckm, int _nParticles,
                        int *_state,
                        double *_probA, double *_probB, double *_rateA, double *_rateB);

void deleteSckmSwarm(SckmSwarm *swarm);
                        
int renormalize(SckmSwarm *swarm);

#endif // SCKM_SWARM_H
