#ifndef __PL_DISCRETE_H__ // guard pl-discrete.h
#define __PL_DISCRETE_H__

#include <R.h>

void cond_discrete_sim_step(const int *nSpecies, const int *nRxns, const int *anStoich,  
                       const double *adHazard, const int *anY, const double *adP, const int *nWhileMax,
                       int *anRxns, int *anX);
void discrete_particle_update(const int *nSpecies, const int *nRxns, const int *anPre, const int *anStoich, 
                              const int *anY, const double *dTau, const *nWhileMax,
                              int *anX, double *adHyper);
void discrete_all_particle_update(const int *nSpecies, const int *nRxns, const int *anPre, const int *anStoich, 
                                  const int *anY, const double *dTau,
                                  const int *nParticles, const int *nWhileMax,
                                  int *anX, double *adHyper);

#endif                   // guard pl-discrete.h
