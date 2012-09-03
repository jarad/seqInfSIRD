
#ifndef GP_H
#define GP_H

typedef struct StochasticChemicalKineticModel {
    unsigned short int s;     // number of species/states
    unsigned short int r;     // number of reactions/transitions
    unsigned short int *Pre;  // r x s matrix of reactants
    unsigned short int *Post; // r x s matrix of products
    short int *Stoich;        // s x r stoichiometry matrix equal to (Post-Pre)'
    double *rates;            // r vector of reaction rate constants
    unsigned int *state;      // s vector of system state
} Sckm;

Sckm *newSckm(const unsigned short int s, const unsigned short int r,
              const unsigned short int *Pre, const unsigned short int *Post,
              const double *rates, const unsigned int *state);


void deleteSckm(Sckm *sckm);

#endif // GP_H
