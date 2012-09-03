#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "Sckm.h"



Sckm *newSckm(const unsigned short int s, const unsigned short int r,
              const unsigned short int *Pre, const unsigned short int *Post,
              const double *rates, const unsigned int *state)
{
    Sckm *sckm;
    sckm = (Sckm *) malloc(sizeof(Sckm));

    sckm->s = s;
    sckm->r = r;

    sckm->Pre = (unsigned short int *) malloc(s*r*sizeof(unsigned short int));
    memcpy(sckm->Pre, Pre, s*r*sizeof(unsigned short int));

    sckm->Post = (unsigned short int *) malloc(s*r*sizeof(unsigned short int));
    memcpy(sckm->Post, Post, s*r*sizeof(unsigned short int));

    sckm->Stoich = (short int *) malloc(s*r*sizeof(short int));
    for (int i=0; i<s; i++) 
    {
        for (int j=0; j<r; j++) 
        {
            sckm->Stoich[i*r+j] = Post[i+j*s]-Pre[i+j*s];
        }
    }

    sckm->rates = (double *) malloc(r*sizeof(double));
    memcpy(sckm->rates, rates, r*sizeof(double));
   
    sckm->state = (unsigned int*) malloc(s*sizeof(unsigned int));
    memcpy(sckm->state, state, s*sizeof(double));

    return(sckm);
}

void deleteSckm(Sckm *sckm)
{
    assert(sckm);
    assert(sckm->Pre);    free(sckm->Pre);
    assert(sckm->Post);   free(sckm->Post);
    assert(sckm->Stoich); free(sckm->Stoich);
    assert(sckm->rates);  free(sckm->rates);
    assert(sckm->state);  free(sckm->state);
    free(sckm);
}


