#ifndef PL_DISCRETE_H // guard pl-discrete.h
#define PL_DISCRETE_H

#include "Sckm.h"
#include "SckmSwarm.h"
#include "SckmParticle.h"

void calc_log_pred_like_R(const int *, const double *, 
                          int *, int *, int *, int *, double *,
                          int *, double *, double *, double *, double *,
                          double *, double *, double *);
double calc_log_pred_like(const int *, double , Sckm *, SckmParticle *, double *);

int cond_discrete_sim_step(Sckm *sckm, const double *adHazard, const int *anY, 
                           const double *adP, int nWhileMax, int *anRxns, int *anX);
int discrete_particle_update(Sckm *sckm, const int *anY, double dTau, int nWhileMax,
                              int *anX, double *adHyper, int *);

void discrete_all_particle_update_R(int *nSpecies, int *nRxns, int *anPre, int *anPost, double *adlMult,
                                  const int *anY, const double *dTau,
                                  int *nParticles, int *nWhileMax,
                                  int *anX, double *adHyper, int *);
int discrete_all_particle_update(Sckm *sckm, 
                                  const int *anY, double dTau,
                                  int nParticles, int nWhileMax,
                                  int *anX, double *adHyper, int *);


void tlpl_R(
           /* Data */
           int *nObs,
           int *anY, 
           double *adTau, 

           /* sckm */
           int *nSpecies, 
           int *nRxns, 
           int *anPre, 
           int *anPost,
           double *adlMult,

           /* Particles */
           int *nParticles, 

           /* Auxiliary */
           int *nResamplingMethod,
           int *nNonuniformity,
           double *dThreshold,
           int *nVerbose,

           /* Outputs */
           int *anX,
           double *adProbA,
           double *adProbB,
           double *adRateA,
           double *adRateB
           );

int tlpl(int nObs, int *anY, double *adTau,
         Sckm *sckm, SckmSwarm **swarm,
         int nResamplingMethod, int nNonuniformity, double dTreshold, int nVerbose);

#endif                   // guard pl-discrete.h
