#ifndef __GILLESPIE_H__ // guard for simulate.h
#define __GILLESPIE_H__ 

void hazard_part_wrap(int *, int *, const int *, const int *, int *);
void hazard_part(int nSpecies, int nRxns, const int *anPre, const int *anX, int *anHp);

void hazard(int nSpecies, int nRxns, const int *anPre, const double *adTheta,   
            const int *anX, double dTau, int *anHp, double *adH); 
void update_species(int nSpecies, int nRxns, const int *anStoich, const int *anRxns,    
                    int *anX);
void tau_leap_one_step(int nSpecies, int nRxns, const int *anStoich, const double *adH,                 
                       int nWhileMax, int *anRxns, int *anX);
void tau_leap_wrap(int *nSpecies, int *nRxns, const int *anStoich, const int *anPre, const double *adTheta,
         const double *adTau, int *nSteps, 
         int *nWhileMax,
         int *anX);
void tau_leap(int nSpecies, int nRxns, const int *anStoich, const int *anPre, const double *adTheta,
         double adTau, int nSteps, int nWhileMax, int *anX); 

int next_to_fire(int nRxns, double *adCuSum);



void gillespie_one_step_wrap(int *nSpecies, int *nRxns, const int *anStoich, const int *anPre, const double *adTheta,
                             double *dT, int *anRxnCount, int *anX);

int gillespie_one_step(int nSpecies, int nRxns, const int *anStoich, const int *anPre, const double *adTheta,
               double dT, int *anRxnCount, int *anX);
void gillespie_wrap(int *nSpecies, int *nRxns, const int *anStoich, const int *anPre, const double *adTheta,
               double *adT, int *nSteps, int *anX);

int gillespie(int nSpecies, int nRxns, const int *anStoich, const int *anPre, const double *adTheta,
               double *adT, int nSteps, int *anX);


#endif                 // guard for simulate.h
