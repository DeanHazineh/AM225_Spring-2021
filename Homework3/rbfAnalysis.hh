#ifndef RBFANALYSIS_HH
#define RBFANALYSIS_HH

#include "rbf.hh"
#include <cstdio>
#include <cmath>
#include <iostream> 
#include<tuple> // for tuple  -- hilbert calculation
#include<vector> // for vector -- hilbert calculation
#include <algorithm> // for sorting functionality
#include "lp_solve.hh"
#include "omp.h"

class rbf_test : rbf{
	public: 	
	/* constructor function-- simply creates Rycroft rbf and conj grad routine indirectly */
		rbf_test(int n_,int type_);
	
	/* initialize point locations */	
		void initPoints(int mode, int hilbertM);
		void demo_LapackdenseSolver();
		void demo_JacobiPCG(int bls);
		double timeSolve_JacobiPCG(int mode);
		void printMatrix(double* mat, int n_);
		void printVector(double* vec, int n_);
		
	/* Count the number of non-zero elements in the dense matrix and the preconditioned matrix */
		int countNonzero(double* mat, int n_);
		int countDenseT();
		int countJacobiP(bool print, int bls);
	
	/* Create initial x,y data points and hold in order sorted by Hilbert curve space */
		void init_Hilbert(int orderM);
		void init_unsortedrandom();
		void print_HilbertData();
		
	private:
		int calculateq(int M, double xval_, double yval_);
		void rot (int N, int &xval, int &yval, int rx, int ry);

};

#endif
