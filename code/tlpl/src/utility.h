
#ifndef __UTILITY_H__ // utility.h guard
#define __UTILITY_H__

int anyNegative(int n, const int *v);                     
void copy(int n, const int *from, int *to);
void runif_vec( int, const double *, const double *, double *);;
void rpois_vec( int, const double *, int *);
void rbeta_vec( int, const double *, const double *, double *);
void rgamma_vec(int, const double *, const double *, double *); 


#endif                // utility.h guard
