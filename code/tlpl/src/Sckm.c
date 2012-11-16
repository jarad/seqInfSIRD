#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "Sckm.h"



Sckm *newSckm(const int s, const int r, int *Pre, int *Post, double *lMult)
{
    Sckm *sckm;
    sckm = (Sckm *) malloc(sizeof(Sckm));

    sckm->s = s;
    sckm->r = r;

    sckm->Pre   = Pre;
    sckm->Post  = Post; 
    sckm->lMult = lMult;

    sckm->Stoich = (int *) malloc(s*r*sizeof(int));
    for (int i=0; i<s; i++) 
    {
        for (int j=0; j<r; j++) 
        {
            sckm->Stoich[i*r+j] = Post[j+i*r]-Pre[j+i*r];
        }
    }

    return(sckm);
}

void deleteSckm(Sckm *sckm)
{
    assert(sckm);
    assert(sckm->Stoich); free(sckm->Stoich);
    free(sckm);
}


