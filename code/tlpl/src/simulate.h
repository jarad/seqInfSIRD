#ifndef __SIMULATE_H__ // guard for simulate.h
#define __SIMULATE_H__ 

void hazard_part(int nSpecies, int nRxns, const int *anPre, const int *anX, int *anHp);
void hazard(int nSpecies, int nRxns, const int *anPre, const double *adTheta,   
            const int *anX, const double *dTau, int *anHp, double *adH); 
void update_species(int nSpecies, int nRxns, const int *anStoich, const int *anRxns,    
                    int *anX);
void sim_one_step(int nSpecies, int nRxns, const int *anStoich, const double *adH,                 
                  int nWhileMax, int *anRxns, int *anX);
void sim(int nSpecies, int nRxns, const int *anStoich, const int *anPre, const double *adTheta,
         const double *adTau, int nSteps, int nWhileMax, int *anX); 


#endif                 // guard for simulate.h
