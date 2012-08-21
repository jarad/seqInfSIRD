#ifndef __PL_DISCRETE_H__ // guard pl-discrete.h
#define __PL_DISCRETE_H__

#include <R.h>

void calculate_log_predictive_likelihood_wrap(int *, int *, const int *, const int *, const int *, const double *,
                                              const double *, double *);
double calculate_log_predictive_likelihood(int , int, const int *, const int *, const int *,
                                           const double *, const double *);

int cond_discrete_sim_step(int nSpecies, int nRxns, const int *anStoich,  
                       const double *adHazard, const int *anY, const double *adP, int nWhileMax,
                       int *anRxns, int *anX);
int discrete_particle_update(int nSpecies, int nRxns, const int *anPre, const int *anStoich, 
                              const int *anY, double dTau, int nWhileMax,
                              int *anX, double *adHyper, int *);

void discrete_all_particle_update_wrap(int *nSpecies, int *nRxns, const int *anPre, const int *anStoich, 
                                  const int *anY, const double *dTau,
                                  int *nParticles, int *nWhileMax,
                                  int *anX, double *adHyper, int *);
int discrete_all_particle_update(int nSpecies, int nRxns, const int *anPre, const int *anStoich, 
                                  const int *anY, double dTau,
                                  int nParticles, int nWhileMax,
                                  int *anX, double *adHyper, int *);

#endif                   // guard pl-discrete.h
