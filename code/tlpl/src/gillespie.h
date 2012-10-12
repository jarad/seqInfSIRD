#include "Sckm.h"

#ifndef GILLESPIE_H // guard for gillespie.h
#define GILLESPIE_H 

void hazard_part_R(int *, int *, int *, int *, const int *, int *);
int hazard_part(Sckm *sckm, const int *anX, int *anHp);

void hazard_R(int *, int *, int *, int *, const double *, const int *, double *, int *, double *);
int hazard(Sckm *sckm, const double *adTheta, const int *anX, double dTau, int *anHp, double *adH); 

void update_species_R(int *, int *, int *, int *, const int *, int *);
int update_species(Sckm *sckm, const int *anRxns, int *anX);

void tau_leap_one_step_R(int *, int *, int *, int *, const double *, int *, int *, int *);
int tau_leap_one_step(Sckm *sckm, const double *adH, int nWhileMax, int *anRxns, int *anX);

void tau_leap_R(int *nSpecies, int *nRxns, int *anPre, int *anPost, 
                   const double *adTheta, const double *adTau, int *nSteps, 
                   int *nWhileMax, int *anRxnCount, int *anX);
int tau_leap(Sckm *sckm, const double *adTheta, const double *adTau, int nSteps, int nWhileMax, int *anRxnCount, int *anX); 

int next_to_fire(int nRxns, double *adCuSum);



void gillespie_one_step_R(int *nSpecies, int *nRxns, int *anPre, int *anPost, 
                             const double *adTheta, double *dT, int *anRxnCount, int *anX);
int gillespie_one_step(Sckm *sckm, const double *adTheta, double dT, int *anRxnCount, int *anX);

void gillespie_R(int *nSpecies, int *nRxns, int *anPre, int *anPost, 
                    const double *adTheta, double *adT, int *nSteps, int *anX);
int gillespie(Sckm *sckm, const double *adTheta, double *adT, int nSteps, int *anX);


#endif                 // guard for gillespie.h
