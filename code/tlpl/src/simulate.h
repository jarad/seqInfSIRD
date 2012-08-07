#ifndef __SIMULATE_H__ // guard for simulate.h
#define __SIMULATE_H__ 

void hazard_part(const int *nSpecies, const int *nRxns, const int *anPre, const int *anX, int *anHp);
void hazard(const int *nSpecies, const int *nRxns, const int *anPre, const double *adTheta,   
            const int *anX, const double *dTau, int *anHp, double *adH); 
void rpois_vec(const int *nRxns, const double *adH, int *anRxns);
void update_species(const int *nSpecies, const int *nRxns, const int *anStoich, const int *anRxns,    
                    int *anX);
void sim_one_step(const int *nSpecies, const int *nRxns, const int *anStoich, const double *adH,                 
                  const int *nWhileMax, int *anRxns, int *anX);
void sim(const int *nSpecies, const int *nRxns, const int *anStoich, const int *anPre, const double *adTheta,
         const double *dTau, const int *nSteps, const int *nWhileMax, int *anX); 


#endif                 // guard for simulate.h
