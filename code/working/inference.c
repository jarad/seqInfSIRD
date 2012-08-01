/*
 *  inference.c
 *  
 *
 *  Created by Jarad Niemi on 2011/2/11
 *  Updated by Jarad Niemi on 2012/6/13
 *
 */

#include <R.h>
#include <Rmath.h>
#include "inference.h"

void hazard(int *anX, int *nN, double *adHazard) {
	adHazard[0] = (double) anX[0]*anX[1]/ *nN;
	adHazard[1] = (double) anX[1];
	adHazard[2] = (double) anX[0];
	adHazard[3] = (double) anX[1];
}

void simulate_one_step(int *anX, double *adTheta, int *nN, int *nRxns, int *anDX) {
	double adHazard[*nRxns];
	hazard(anX, nN, adHazard);
	
	int i;
	for (i=0; i<*nRxns; i++) {
		anDX[i] = rpois(adTheta[i]*adHazard[i]);
	}
	
	// Restrict to maintain non-negative states
	anDX[0] = anDX[0] < anX[0] ? anDX[0] : anX[0];
	if (anDX[1]+anDX[3] > anX[1]+anDX[0]) { // split I-> rxns between I->R and I->D
		anDX[1] = rbinom(anX[1]+anDX[0], adHazard[1]/(adHazard[1]+adHazard[3]));
		anDX[3] = anX[1]+anDX[0]-anDX[1];
	}
	anDX[2] = anDX[2] < anX[0]-anDX[0] ? anDX[2] : anX[0]-anDX[0];
	
	// Update states
	anX[0] += -anDX[0]        -anDX[2]        ;
	anX[1] +=  anDX[0]-anDX[1]        -anDX[3];
	anX[2] +=          anDX[1]+anDX[2]        ;
	anX[3] +=                         +anDX[3];
	
	if (anX[0]<0 || anX[1]<0 || anX[2]<0 || anX[3]<0) {
		Rprintf("Negative values for states!!!!\n");
	}
}

void simulate(int *anX, double *adTheta, int *nN, int *nSteps, int *nStates, int *nReps, int *nRxns, int *anDX) {
	int i, j, k, nOffset=0;
	
	GetRNGstate();
	for (i=0; i<*nReps; i++) {
		for (j=0; j<*nSteps; j++) {
			nOffset += *nStates;
			for (k=0; k<*nStates; k++) anX[nOffset+k] = anX[nOffset+k-(*nStates)]; // copy states over
			simulate_one_step(&anX[nOffset], adTheta, &nN[i], nRxns, &anDX[nOffset]);
		}
		nOffset += *nStates; // get to new starting state
	}
	PutRNGstate();
}


void inference_one_step(int *anX0, int *anDX, int *nN, double *adProp, double *adHyper, int *nRxns, int *bSample) {
	double adHazard[*nRxns];
	hazard(anX0, nN, adHazard);
	
	int i, offset;
	for (i=0; i<*nRxns; i++) {
		if (*bSample) {
			adHyper[i] += rbinom( anDX[i],adProp[i]); 
		} else {
			adHyper[i] += adProp[i] * anDX[i];
		}
		adHyper[i+ *nRxns] += adProp[i] * adHazard[i]; 
	}
}


void inference(int *anX0, int *anDX, int *nN, double *adProp, double *adHyper, int *nRxns, int *nStates, int *nReps) {
	int i, j, nOffset=0;
	for (i=0; i<*nReps; i++) {
		// doesn't do anything yet
	}
}


void one_step(int *anX, double *adHyper, double *adTheta, double *adProp, int *nN, 
              int *nReps, int *nStates, int *nRxns, int *bSample, int *anDX) {
	int i, j, nOffsetStates=0, nOffsetRxns=0, anX0[*nStates];
	
	GetRNGstate();
	for (i=0; i< *nReps; i++) {
		for (j=0; j<*nStates; j++) anX0[j] = anX[nOffsetStates+j];
		simulate_one_step(&anX[nOffsetStates], &adTheta[nOffsetRxns], &nN[i], nRxns, &anDX[nOffsetRxns]);
		inference_one_step(anX0, &anDX[nOffsetRxns], &nN[i], &adProp[nOffsetRxns], &adHyper[nOffsetRxns*2], nRxns, bSample);
		
		nOffsetStates += *nStates;
		nOffsetRxns   += *nRxns;
	}
	PutRNGstate();
}

