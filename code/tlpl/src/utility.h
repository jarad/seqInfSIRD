
#ifndef __UTILITY_H__ // utility.h guard
#define __UTILITY_H__

int anyNegative(int n, int *v);                     
void copy(int n, int *from, int *to);
void rpois_vec(const int *nRxns, const double *adH, int *anRxns);
void rbeta_vec( const int *, const double *, const double *, double *);
void rgamma_vec(const int *, const double *, const double *, double *); 


#endif                // utility.h guard
