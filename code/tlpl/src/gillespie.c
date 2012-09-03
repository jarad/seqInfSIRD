/* 
* Functions for discrete-time compartment models
*/

#include <R.h>
#include <Rmath.h>
#include "utility.h"
#include "gillespie.h"
#include "Sckm.h"


/* Calculates the part of the hazard other than the fixed parameter */
void hazard_part_wrap(int *nSpecies, int *nRxns, int *anPre, int *anPost, const int *anX, int *anHazardPart)
{
    Sckm *sckm = newSckm(*nSpecies, *nRxns, anPre, anPost);
    hazard_part(sckm, anX, anHazardPart);
    deleteSckm(sckm);
}

int hazard_part(Sckm *sckm, const int *anX, int *anHazardPart)
{
    int nSpecies=sckm->s;
    for (int i=0; i<sckm->r; i++) 
    {
        anHazardPart[i] = 1;
        for (int j=0; j<nSpecies; j++) 
        {
            anHazardPart[i] *= choose(anX[j], sckm->Pre[i* nSpecies+j]); 
        }    
    }  
    return 0;     
}

/* Calculates the hazard for the next reaction */
void hazard_wrap(int *nSpecies, int *nRxns, int *anPre, int *anPost,
            const double *adTheta,   
            const int *anX, 
            double *dTau,  
            int *anHazardPart, double *adHazard) 
{
    Sckm *sckm = newSckm(*nSpecies, *nRxns, anPre, anPost);
    hazard(sckm, adTheta, anX, *dTau, anHazardPart, adHazard);
    deleteSckm(sckm);
}


int hazard(Sckm *sckm, const double *adTheta, const int *anX, double dTau,  
            int *anHazardPart, double *adHazard)                   // return: hazard part and hazard
{
    hazard_part(sckm, anX, anHazardPart);
    for (int i=0; i<sckm->r; i++) 
    {
        adHazard[i] = adTheta[i]*anHazardPart[i]*dTau;
    }
    return 0;
}


/* Updates the species according to the stoichiometry */
void update_species_wrap(int *nSpecies, int *nRxns, int *anPre, int *anPost,   
                    const int *anRxnCount, int *anX)
{
    Sckm *sckm = newSckm(*nSpecies, *nRxns, anPre, anPost);
    update_species(sckm, anRxnCount, anX);
    deleteSckm(sckm);
}

int update_species(Sckm *sckm, const int *anRxnCount, int *anX)              // return: updated species
{
    int i,j;
    for (i=0; i<sckm->s; i++) 
    {
        for (j=0; j<sckm->r; j++) 
        {
            anX[i] += sckm->Stoich[sckm->s*j+i]*anRxnCount[j];
        }    
    } 
    return 0;
}




// --------------------------------------------------------------------------------------
// Discrete time
// --------------------------------------------------------------------------------------


/* Forward simulate ahead one time-step */
void tau_leap_one_step_wrap(int *nSpecies, int *nRxns, int *anPre, int *anPost, 
                  const double *adHazard,                 
                  int *nWhileMax,                          
                  int *anRxnCount, int *anX)
{
    Sckm *sckm = newSckm(*nSpecies, *nRxns, anPre, anPost);
    tau_leap_one_step(sckm, adHazard, *nWhileMax, anRxnCount, anX);
    deleteSckm(sckm);
}

int tau_leap_one_step(Sckm *sckm, 
                  const double *adHazard,                 
                  int nWhileMax,                          
                  int *anRxnCount, int *anX)                                 // return: updated species
{
    int nSpecies=sckm->s;
    int i, whileCount=0, anTempX[nSpecies];
    while (1) 
    {
        memcpy(anTempX, anX, nSpecies*sizeof(int));

        // Get number of reactions
        int i;
        GetRNGstate();
        for (i=0; i<sckm->r; i++) anRxnCount[i] = rpois(adHazard[i]);
        PutRNGstate();

        update_species(sckm, anRxnCount, anTempX);

        if (!anyNegative(nSpecies, anTempX)) 
        {
            memcpy(anX, anTempX, nSpecies*sizeof(int));
            return 0;
        }

        // Limit how long the simulation tries to find a non-negative update
        whileCount++;
        if (whileCount>nWhileMax)  
        {
            error("C:tau_leap_one_step: Too many unsuccessful simulation iterations.");
            return 1;
        }
    }
} 



void tau_leap_wrap(int *nSpecies, int *nRxns, int *anPre, int *anPost, const double *adTheta,
         const double *adTau, int *nSteps, 
         int *nWhileMax,
         int *anX)
{
    Sckm *sckm = newSckm(*nSpecies, *nRxns, anPre, anPost);
    tau_leap(sckm, adTheta, adTau, *nSteps, *nWhileMax, anX);
    deleteSckm(sckm);
}

int tau_leap(Sckm *sckm, const double *adTheta,
         const double *adTau, int nSteps, 
         int nWhileMax,
         int *anX)
{
    int nRxns = sckm->r, nSpecies = sckm->s;
    int i, *ipLast, *ipCurrent, anRxnCount[nRxns], anHazardPart[nRxns];
    ipLast     = anX;             // Points to last state
    ipCurrent  = ipLast+nSpecies; // Points to current state

    double adHazard[nRxns];
    for (i=0; i<nSteps;i++)
    {
        memcpy(ipCurrent, ipLast, nSpecies*sizeof(int));
        hazard(sckm, adTheta, ipCurrent, adTau[i], anHazardPart, adHazard);
        tau_leap_one_step(sckm, adHazard, nWhileMax, anRxnCount, ipCurrent);
        ipLast += nSpecies; ipCurrent += nSpecies;
    }
    return 0;
}




// --------------------------------------------------------------------------------------
// Continuous time
// --------------------------------------------------------------------------------------


int next_to_fire(int nRxns, double *adCuSum) {
    int i, next=0;
    double dUniform, dSum=adCuSum[nRxns-1]; 
    
    GetRNGstate();
    dUniform = runif(0,1);
    PutRNGstate();

    for (i=0; i<nRxns; i++) {
        adCuSum[i] /= dSum;
        if (dUniform < adCuSum[i]) return next;
        next++;
    }   
    // print error if this fails
}


void gillespie_one_step_wrap(int *nSpecies, int *nRxns, int *anPre, int *anPost, 
                             const double *adTheta, double *dT,  int *anRxnCount, int *anX)
{
    Sckm *sckm = newSckm(*nSpecies, *nRxns, anPre, anPost);
    gillespie_one_step(sckm, adTheta, *dT, anRxnCount, anX);
    deleteSckm(sckm);
}


int gillespie_one_step(Sckm *sckm, const double *adTheta, double dT, int *anRxnCount, int *anX)
{
    int nSpecies = sckm->s, nRxns = sckm->r;
    int i, nRxnID, anX0[nSpecies], anHazardPart[nRxns];
    double dCurrentTime=0, adHazard[nRxns], adCuSum;
    while (1) {
        memcpy(anX0, anX, nSpecies*sizeof(int));
        hazard(sckm, adTheta, anX, 1, anHazardPart, adHazard);
        
        // Calculate cumulative hazard
        for (i=1; i<nRxns; i++) adHazard[i] += adHazard[i-1];
        if (adHazard[nRxns-1] < 0.0001) return 0;                 // make this a function of dT?

        dCurrentTime += rexp(1/adHazard[nRxns-1]);
        if (dCurrentTime > dT) return 0;                          // stopping condition
       
        nRxnID = next_to_fire(nRxns, adHazard);
        anRxnCount[nRxnID]++;

        for (i=0; i<nSpecies; i++) anX[i] += sckm->Stoich[nSpecies * nRxnID + i]; 
    }
    return 0;               
}



void gillespie_wrap(int *nSpecies, int *nRxns, int *anPre, int *anPost, const double *adTheta,
               double *adT, int *nSteps, int *anX)
{
    Sckm *sckm = newSckm(*nSpecies, *nRxns, anPre, anPost);
    gillespie(sckm, adTheta, adT, *nSteps, anX);
    deleteSckm(sckm);
}




int gillespie(Sckm *sckm, const double *adTheta, double *adT, int nSteps, int *anX)
{
    int i, nSO=0, anRxnCount[sckm->r], anHazardPart[sckm->r];
    double adHazard[sckm->r];
    for (i=0; i<nSteps;i++)
    {
        memcpy(&anX[nSO+sckm->s], &anX[nSO], sckm->s*sizeof(int));
        nSO += sckm->s;

        gillespie_one_step(sckm, adTheta, adT[i], anRxnCount,  &anX[nSO]);
    }
    return 0;
}


