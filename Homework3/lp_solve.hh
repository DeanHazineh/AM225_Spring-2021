#ifndef LP_SOLVE_HH
#define LP_SOLVE_HH

void lapack_fail();
void solve_banded_matrix(int n,int kl,int ku,double *A,double *x);


#endif
