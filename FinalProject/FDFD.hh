// Dean Hazineh
// AP225 Specific
// Creates the Yee's 3D FD Grid for defining central difference quantities
// which is reformulated into the linear problem Ae=b;

// update--actually I don't think it does implement the Yee cell concept
//--think the implementation is wrong here! Not sure at all

#ifndef FDFD_HH
#define FDFD_HH
#include "blas.h"
#include "cblas.h"
#include <cstdio>
#include <cmath>
#include <complex>// std::complex
typedef std::complex<double> dcomplex;
#include <stdexcept>

class FDFD_grid{
	public: 
		// Test parameter for epsilon object
		double sphereRad = 3e-6;
	
		// Define the frequency for the problem 
		const double omega;
		// sigma double used for SC-PML
		const double m; //Polynomial PML Order
		double sigma;
		
		// Specify the number of points for grid region and PML ghost nodes
		const int nx;
		const int ny; 
		const int nz;
		const int nPML;
		const int fieldOriginIdx;
		const int zsliceGridSize;
		// Specify the uniform grid spacing and the PML discretization
		const double dx;
		const double dy;
		const double dz;
		const double dPML;
		// Hold the resulting lenghth in each grid direction and PML on one side
		const double Lx;
		const double Ly; 
		const double Lz;
		const double LPML;
		/* gridLen is length of solution/source row vector with n row mapped to (w,i,j,k) where w=x,y,z (length includes the ghost nodes)*/
		const int gridLen; 
		
		FDFD_grid(double omega_, double dx_, double dy_, double dz_, double dPML_, int nx_, int ny_, int nz_, int nPML_);
		~FDFD_grid();
		void fillSpaceMaterial();
		void viewEps();
		void defineSourceTerm();
		void viewJ();
		void AssembleAMatrixDoubleDemo();
		void AssembleAMatrix();
		void BICGSolve(); 
		void viewESol();
		
	private:
		const double epsilonnaught = 8.854187817e-12;
		const double munaught = 1.25663706212e-6;
		// hold the electric field solution --row vector
		dcomplex* e;
		// hold the effective current source density-- row vector
		dcomplex* j;
		// hold the magnetic and electric permeability diagonal elements -- Operators
		//Dmu and Deps are just diag_mu*I and diag_eps*I 
		dcomplex* diag_mu;
		dcomplex* diag_eps;

		// hold the intermediate curlE operator and curlH operator
		// Note defining dcomplex structure for matrix gridLen*gridLen causes stack overflow
		// for any reasonable grid size-- plus the curl itself is double real 
		double* Ch_real;
		double* Ce_real;
		double* Ch_imag;
		double* Ce_imag;
		double* Areal;
		double* Aimag;
		double* interm;
		double* interm2;
		dcomplex* A; 
			
		void output3D(dcomplex* dataVector, const char* fileOutput);	
		void displayTransferMatrix(double* matrix, const char* fileOutput);	
		void BlasMultiply(int N, double* A, double* B, double* C, double alpha, double beta);
		void matrixAdd(int N, double* A, double* B, double* C, double alpha);
//		void BLAS_cgemv(const char trans, int m, int n, dcomplex* alpha, dcomplex *a, int lda, dcomplex* x, int incx,
//				dcomplex* beta, dcomplex *y, int incy);
		void complexConjugateVector(int N, dcomplex* a, dcomplex* b);
		
		const char transn = 'n';
		const char transc = 'c';
		const char transt = 't';
		//dcomplex* identityMatrix; 

};


#endif












