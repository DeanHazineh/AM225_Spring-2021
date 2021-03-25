#include "lp_solve.hh"

#include <cstdio>
#include <cstdlib>

struct doublecomplex
{   double re;
   double im;
}; 

// Tell the compiler about the existence of the required LAPACK functions
extern "C" {
    // the zgbsv solver which is the complex band matrix solver
    // http://www.netlib.org/lapack/explore-html/d9/dbb/group__complex16_g_bsolve_ga908abc0aad64131b9a32edb08510eb00.html
    void zgbsv_(int *n,int *kl,int *ku,int *nrhs,doublecomplex *ab,int *ldab,int *ipiv,
            doublecomplex *b,int *ldb,int *info);
}

/** Prints a generic error message in cases when LAPACK reports an error. */
void lapack_fail() {
    fputs("LAPACK routine failed\n",stderr);
    exit(1);
}


/** Solves the matrix system Ax=b for the case when A is banded.
 * \param[in] n the dimension of the matrix.
 * \param[in] (kl,lu) the number of subdiagonals and superdiagonals.
 * \param[in] A the matrix terms; see DGBSV documentation for full details.
 * \param[in] x the source vector, which will be replaced by the solution
 *              when the routine exits. */
void solve_banded_matrix(int n,int kl,int ku,double *A,double *x) {

    // Create the temporary memory that LAPACK needs
    int info,nrhs=1,*ipiv=new int[n],ldab=2*kl+ku+1;

    // Make call to LAPACK routine
    zgbsv_(&n,&kl,&ku,&nrhs,A,&ldab,ipiv,x,&n,&info);
    if(info!=0) {
        fputs("LAPACK routine failed\n",stderr);
        exit(1);
    }

    // Remove temporary memory
    delete [] ipiv;
}

