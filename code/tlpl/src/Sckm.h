/* 
Sckm class definition. Defines a STOCHASTIC CHEMICAL KINETIC MODEL consisting of
  - s: number of species
  - r: number of reactions
  - pre: r x s matrix indicating species consumed in a reaction and
                      determine the reaction order
  - post: r x s matrix indicating species produced in a reaction
  - stoich: s x r matrix equal to the transpose of (Post-Pre)
  - rate_constant: r vector of rate constants
  - state: s vector for current system state 
*/

#ifndef SCKM_H
#define SCKM_H


class Sckm
{
  public:
    Sckm();
    Sckm(unsigned int s, unsigned int r);
    Sckm(unsigned int s, unsigned int r, 
         unsigned int *pre, unsigned int *post, int *stoich, 
         double *rates, unsigned int *state);
    virtual ~Sckm();

    unsigned int getNumSpecies() { return s; }
    unsigned int getNumRxns()    { return r; }


  protected:
    unsigned int s,r, *pre, *post, *state;
    int *stoich;
    double *rates;
    bool allocatedPre, 
         allocatedPost,
         allocatedStoich,
         allocatedRates,
         allocatedState;
};

#endif
