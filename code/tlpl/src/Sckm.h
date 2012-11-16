
#ifndef SCKM_H
#define SCKM_H

typedef struct StochasticChemicalKineticModel {
    int s;         // number of species/states
    int r;         // number of reactions/transitions
    int *Pre;      // r x s matrix of reactants
    int *Post;     // r x s matrix of products
    int *Stoich;   // s x r stoichiometry matrix equal to (Post-Pre)'
    double *lMult; // r vector of constants in reaction rate calculation
} Sckm;

Sckm *newSckm(const int s, const int r, int *Pre, int *Post, double *lMult);

void deleteSckm(Sckm *sckm);

#endif // SCKM_H
