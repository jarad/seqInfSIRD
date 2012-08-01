/*
 *  inference.h
 *  
 *
 *  Created by Jarad Niemi on 2011/2/11
 *  Updated by Jarad Niemi on 2012/6/13
 *
 */

void hazard(int *anX, int *nN, double *adHazard);
void simulate_one_step(int *anX, double *adTheta, int *nN, int *nRxns, int *anDX);
void inference_one_step(int *anX0, int *anDX, int *nN, double *adProp, double *adHyper, int *nRxns, int *bSample);
void inference(int *anX0, int *anDX, int *nN, double *adProp, double *adHyper, int *nRxns, int *nStates, int *nReps);
void one_step(int *anX, double *adHyper, double *adTheta, double *adProp, int *nN, int *nReps, int *nStates, int *nRxns, int *bSample, int *anDX);

