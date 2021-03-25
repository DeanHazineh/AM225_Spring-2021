#ifndef STRASSENFUNCTIONS_HH
#define STRASSENFUNCTIONS_HH
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <omp.h>
#include <algorithm>    // std::copy
#include "blas.h"

// Generate a pseudo random double between 0 and upperRange
	inline double doubleRandGen(double upperRange){
		return (double)rand()/RAND_MAX*upperRange;
	}

// Print matrix to the terminal
	void printMat(int N, double* mat);

// Generate 2D matrix (NxN) with random elements
	double* genRandomMat(int N);

// Do standard 2D Matrix Multiplication 
// with three for loops and pulling elements from heap arrays
	double* mat_mul(int N, double* A, double* B);
// Try the standard 2D Matrix Multiply but 
// with storing segments into fast memory for arithmatic
	double* mat_mul_menManage(int N, double* A, double* B);

// STRASSEN ROUTINE
// Split matrix into Submatrices as required for Strassen Algorithm
	void splitStrassen(int N, double* in, double* zerozero, double* zeroone, double* onezero, double* oneone);

// Define Matrix Addition A + alpha B
	void matAdd(int N, double* out, double* A, double* B, double alpha);
	void matAdd_memManage(int N, double* out, double* A, double* B, double alpha);

// Recombine submatrices into complete matrix as required for strassen (inverse operation of splitStrassen)
	void recombineStrassen(int N, double* out, double* Czerozero, double* Czeroone, double* Conezero, double*Coneone);

// conduct matrix mulitplication via the Strassen Algorithm
	void writeStrassen(int N, double* out, double* A, double* B);
	double* writeStrassen2(int N, double* A, double* B);

// BLAS 
	double* BlasMultiply(int N, double* A, double* B);



	
#endif
