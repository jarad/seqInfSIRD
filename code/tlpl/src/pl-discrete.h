#ifndef PL_DISCRETE_H // guard pl-discrete.h
#define PL_DISCRETE_H

#include "Sckm.h"

void calc_log_pred_like_wrap(int *, int *, int *, int *, 
                             const int *, const int *, const double *, const double *, double *);
double calc_log_pred_like(Sckm *sckm, const int *, const int *, const double *, const double *);

int cond_discrete_sim_step(Sckm *sckm, const double *adHazard, const int *anY, 
                           const double *adP, int nWhileMax, int *anRxns, int *anX);
int discrete_particle_update(Sckm *sckm, const int *anY, double dTau, int nWhileMax,
                              int *anX, double *adHyper, int *);

void discrete_all_particle_update_wrap(int *nSpecies, int *nRxns, int *anPre, int *anPost, 
                                  const int *anY, const double *dTau,
                                  int *nParticles, int *nWhileMax,
                                  int *anX, double *adHyper, int *);
int discrete_all_particle_update(Sckm *sckm, 
                                  const int *anY, double dTau,
                                  int nParticles, int nWhileMax,
                                  int *anX, double *adHyper, int *);

#endif                   // guard pl-discrete.h
