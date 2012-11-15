#ifndef SCKM_SWARM_H
#define SCKM_SWARM_H

#include "Sckm.h"
#include "SckmParticle.h"

typedef struct StochasticChemicalKineticModelSwarm {
    int nParticles, nStates, nRxns, logWeights, normalizedWeights;
    double *dWeights;
    SckmParticle *pParticle;
} SckmSwarm;

SckmSwarm *newSckmSwarm(Sckm *sckm, int _nParticles,
                        int *_state,
                        double *_probA, double *_probB, double *_rateA, double *_rateB);

void deleteSckmSwarm(SckmSwarm *swarm);
                        
int renormalize(SckmSwarm *swarm);

#endif // SCKM_SWARM_H
