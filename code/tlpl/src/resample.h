
int compare_doubles(const void *, const void *);

void is_unsorted_wrap(int *, const double *, int*);
int is_unsorted(int , const double *);

void cumulative_sum_wrap(int *, double *);
int cumulative_sum(int , double *);

void rep2id_wrap(int *, int *, int *);
int rep2id(int *, int , int *);

void inverse_cdf_weights_wrap(int *, double *, int *, double *, int *);
int inverse_cdf_weights(int , double *, int , double *, int *);


void effective_sample_size_wrap(int *, double *, double *);
double effective_sample_size(int , double *);

void coefficient_of_variation_wrap(int *, double *, double *);
double coefficient_of_variation(int , double *);

void entropy_wrap(int *, double *, double *);
double entropy(int , double *);

void renormalize_wrap(int *, int *, double *);
int renormalize(int, int, double *);

void multinomial_resample_wrap(int *, double *, int *, int *);
int multinomial_resample(int, double *, int, int *);

void stratified_resample_wrap(int *, double *, int *, int *);
int stratified_resample(int , double *, int , int *);

void systematic_resample_wrap(int *, double *, int *, int *);
int systematic_resample(int , double *, int , int *);

void residual_resample_wrap(int *, double *, int *, int *, int *);
int residual_resample(int , double *, int , int *, int );


