#ifndef PL_DISCRETE_H // guard pl-discrete.h
#define PL_DISCRETE_H

#include "Sckm.h"
#include "SckmParticle.h"

void calc_log_pred_like_R(const int *, const double *, 
                          int *, int *, int *, int *, 
                          int *, double *, double *, double *, double *,
                          double *, double *, double *);
double calc_log_pred_like(const int *, double , Sckm *, SckmParticle *);

int cond_discrete_sim_step(Sckm *sckm, const double *adHazard, const int *anY, 
                           const double *adP, int nWhileMax, int *anRxns, int *anX);
int discrete_particle_update(Sckm *sckm, const int *anY, double dTau, int nWhileMax,
                              int *anX, double *adHyper, int *);

void discrete_all_particle_update_R(int *nSpecies, int *nRxns, int *anPre, int *anPost, 
                                  const int *anY, const double *dTau,
                                  int *nParticles, int *nWhileMax,
                                  int *anX, double *adHyper, int *);
int discrete_all_particle_update(Sckm *sckm, 
                                  const int *anY, double dTau,
                                  int nParticles, int nWhileMax,
                                  int *anX, double *adHyper, int *);

#endif                   // guard pl-discrete.h
