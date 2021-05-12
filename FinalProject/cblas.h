#ifndef CBLAS_H
#define CBLAS_H

#include <complex>   // std::complex
typedef std::complex<double> dcomplex;
#define F77_CALL(x)  x ## _
#define F77_NAME(x)  F77_CALL(x)

#ifdef  __cplusplus
extern "C" {
#endif

	extern int zgbsv_(const int *n, const int *kl, const int *ku, const int *nrhs,
		            dcomplex *ab, const int *ldab, int *ipiv, dcomplex *b,
		            const int *ldb, int *info);

	extern void cgemv_(const char *trans, const int *m, const int *n,
		                           dcomplex *alpha, dcomplex *a,
		                           const int *lda, dcomplex *x, const int *incx,
		                           dcomplex *beta, dcomplex *y, const int *incy);


	extern void cgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k,
		                           dcomplex *alpha, dcomplex *a, const int *lda,
		                           dcomplex *b, const int *ldb,
		                           dcomplex *beta, dcomplex *c, const int *ldc);      
		                           

	extern void ccopy_(const int *n, const dcomplex *dx, const int *incx,
		                           dcomplex *dy, const int *incy);      
		   
	extern dcomplex cdotu_(const int *n, dcomplex* cx, const int* incx, dcomplex* cy, const int* incy);
	extern dcomplex cdotc_(const int *n, dcomplex* cx, const int* incx, dcomplex* cy, const int* incy);
	
	extern void caxpy_(const int *n, dcomplex* ca, dcomplex* cx, const int* incx, dcomplex* cy, const int* incy);
	
		
#ifdef  __cplusplus
}
#endif

#endif /* BLAS_H */



