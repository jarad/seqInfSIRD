/* 
* Functions for discrete-time compartment models
*/

#include <R.h>
#include <Rmath.h>
#include "utility.h"
#include "gillespie.h"

/* Calculates the part of the hazard other than the fixed parameter */
void hazard_part_wrap(int *nSpecies, int *nRxns, const int *anPre, const int *anX, int *anHazardPart)
{
    hazard_part(*nSpecies, *nRxns, anPre, anX, anHazardPart);
}

void hazard_part(int nSpecies, int nRxns, const int *anPre, // system specific arguments 
                 const int *anX,                                          // current system state
                 int *anHazardPart)                                       // return: (partial) system hazard
{
    int i, j;
    for (i=0; i<nRxns; i++) 
    {
        anHazardPart[i] = 1;
        for (j=0; j<nSpecies; j++) 
        {
            anHazardPart[i] *= choose(anX[j], anPre[i* nSpecies+j]); 
        }    
    }       
}

/* Calculates the hazard for the next reaction */
void hazard_wrap(int *nSpecies, int *nRxns, const int *anPre, const double *adTheta,   
            const int *anX, 
            double *dTau,  
            int *anHazardPart, double *adHazard) 
{
    hazard(*nSpecies, *nRxns, anPre, adTheta, anX, *dTau, anHazardPart, adHazard);
}


void hazard(int nSpecies, int nRxns, const int *anPre, const double *adTheta,   
            const int *anX, 
            double dTau,  
            int *anHazardPart, double *adHazard)                   // return: hazard
{
    hazard_part(nSpecies, nRxns, anPre, anX, anHazardPart);
    int i;
    for (i=0; i<nRxns; i++) 
    {
        adHazard[i] = adTheta[i]*anHazardPart[i]*dTau;
    }
//    Rprintf("\n");
}


/* Updates the species according to the stoichiometry */
void update_species_wrap(int *nSpecies, int *nRxns,     
                    const int *anStoich, const int *anRxnCount,    
                    int *anX)
{
    update_species(*nSpecies, *nRxns, anStoich, anRxnCount, anX);
}

void update_species(int nSpecies, int nRxns,     
                    const int *anStoich, const int *anRxnCount,    
                    int *anX)                               // return: updated species
{
    int i,j;
    for (i=0; i<nSpecies; i++) 
    {
        for (j=0; j<nRxns; j++) 
        {
            anX[i] += anStoich[nSpecies*j+i]*anRxnCount[j];
        }    
    } 
}




// --------------------------------------------------------------------------------------
// Discrete time
// --------------------------------------------------------------------------------------


/* Forward simulate ahead one time-step */
void tau_leap_one_step_wrap(int *nSpecies, int *nRxns, const int *anStoich, 
                  const double *adHazard,                 
                  int *nWhileMax,                          
                  int *anRxnCount, int *anX)
{
    tau_leap_one_step(*nSpecies,*nRxns, anStoich, adHazard, *nWhileMax, anRxnCount, anX);
}

void tau_leap_one_step(int nSpecies, int nRxns, const int *anStoich, 
                  const double *adHazard,                 
                  int nWhileMax,                          
                  int *anRxnCount, int *anX)                                 // return: updated species
{
    int  whileCount=0, anTempX[nSpecies];
    while (1) 
    {
        // Copy current state for temporary use
        memcpy(anTempX, anX, nSpecies*sizeof(int));
        //copy(nSpecies, anX, anTempX);

        // Get number of reactions
        rpois_vec(nRxns, adHazard, anRxnCount);

        // Update species
        update_species(nSpecies, nRxns, anStoich, anRxnCount, anTempX);

        // Test if update has any negative species
        if (!anyNegative(nSpecies, anTempX)) 
        {
            // Copy temporary state into current state
            memcpy(anX, anTempX, nSpecies*sizeof(int));
            //copy(nSpecies, anTempX, anX);
            break;
        }

        // Limit how long the simulation tries to find a non-negative update
        whileCount++;
        if (whileCount>nWhileMax) error("C:tau_leap_one_step: Too many unsuccessful simulation iterations.");
    }
} 



void tau_leap_wrap(int *nSpecies, int *nRxns, const int *anStoich, const int *anPre, const double *adTheta,
         const double *adTau, int *nSteps, 
         int *nWhileMax,
         int *anX)
{
    tau_leap(*nSpecies, *nRxns, anStoich, anPre, adTheta, *adTau, *nSteps, *nWhileMax, anX);
}

void tau_leap(int nSpecies, int nRxns, const int *anStoich, const int *anPre, const double *adTheta,
         double adTau, int nSteps, 
         int nWhileMax,
         int *anX)
{
    int i, nSO=0, anRxnCount[nRxns], anHazardPart[nRxns];
    double adHazard[nRxns];
    for (i=0; i<nSteps;i++)
    {
        // Copy state for current step
        memcpy(&anX[nSO+nSpecies], &anX[nSO], nSpecies*sizeof(int));
        //copy(nSpecies, &anX[nSO], &anX[nSO+nSpecies]);
        nSO += nSpecies;

        // Calculate hazard based on current state
        hazard(nSpecies, nRxns, anPre, adTheta, &anX[nSO], adTau, anHazardPart, adHazard);

        // Forward simulate the system
        tau_leap_one_step(nSpecies, nRxns, anStoich, adHazard, nWhileMax, anRxnCount, &anX[nSO]);
    }
}




// --------------------------------------------------------------------------------------
// Continuous time
// --------------------------------------------------------------------------------------


int next_to_fire(int nRxns, double *adCuSum) {
    int i, next=0;
    double dUniform=runif(0,1), dSum=adCuSum[nRxns-1]; 
    for (i=0; i<nRxns; i++) {
        adCuSum[i] /= dSum;
        if (dUniform < adCuSum[i]) return next;
        next++;
    }   
    // print error if this fails
}


void gillespie_one_step_wrap(int *nSpecies, int *nRxns, const int *anStoich, const int *anPre, const double *adTheta,
                             double *dT,  int *anRxnCount, int *anX)
{
    gillespie_one_step(*nSpecies, *nRxns, anStoich, anPre, adTheta, *dT, anRxnCount, anX);
}


int gillespie_one_step(int nSpecies, int nRxns, const int *anStoich, const int *anPre, const double *adTheta,
               double dT, int *anRxnCount, int *anX)
{
    int i, nRxnID, anX0[nSpecies], anHazardPart[nRxns];
    double dCurrentTime=0, adHazard[nRxns], adCuSum;
    while (1) {
        memcpy(anX0, anX, nSpecies*sizeof(int));
        //copy(nSpecies, anX, anX0);
        hazard(nSpecies, nRxns, anPre, adTheta, anX, 1, anHazardPart, adHazard);
        
        // Calculate cumulative hazard
        for (i=1; i<nRxns; i++) adHazard[i] += adHazard[i-1];
        if (adHazard[nRxns-1] < 0.0001) return 0;                 // make this a function of dT?

        dCurrentTime += rexp(1/adHazard[nRxns-1]);
        if (dCurrentTime > dT) return 0;                          // stopping condition
       
        nRxnID = next_to_fire(nRxns, adHazard);
        anRxnCount[nRxnID]++;

        for (i=0; i<nSpecies; i++) anX[i] += anStoich[nSpecies * nRxnID + i]; 
    }
    return 0;               
}



void gillespie_wrap(int *nSpecies, int *nRxns, const int *anStoich, const int *anPre, const double *adTheta,
               double *adT, int *nSteps, int *anX)
{
    gillespie(*nSpecies, *nRxns, anStoich, anPre, adTheta, adT, *nSteps, anX);
}




int gillespie(int nSpecies, int nRxns, const int *anStoich, const int *anPre, const double *adTheta,
               double *adT, int nSteps, int *anX)
{
    int i, nSO=0, anRxnCount[nRxns], anHazardPart[nRxns];
    double adHazard[nRxns];
    for (i=0; i<nSteps;i++)
    {
        // Copy state for current step
        memcpy(&anX[nSO+nSpecies], &anX[nSO], nSpecies*sizeof(int));
        //copy(nSpecies, &anX[nSO], &anX[nSO+nSpecies]);
        nSO += nSpecies;

        // Forward simulate the system
        gillespie_one_step(nSpecies, nRxns, anStoich, anPre, adTheta, adT[i], anRxnCount,  &anX[nSO]);
    }
    return 0;
}


