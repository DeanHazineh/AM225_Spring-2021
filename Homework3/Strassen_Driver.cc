// Code for Problem 1 on Strassen's Algorithm
#include "StrassenFunctions.hh"
#include <cstdio>


void testMatrixMultiply();
void checkBLASStrassenAgainstStandard(int N, bool Print);
void timeStandardMultiply();
void timeStrassenMultiply();
void timeBlas();

int main(){	

	/* Test two versions of standard matrix mulitply and print
		-> In doing this, one observes a huge effect in reading from heap to fast memory
	 	-> vs operating on the heap array
	*/ 
		testMatrixMultiply();
		
	
	/* Print out the standard, Strassen, and Blas matrix mulitply result to validate accuracy
		-> check two different values since the Strassen is recursive with a standard mat mul below certain N
	*/
		checkBLASStrassenAgainstStandard(6, true);
		checkBLASStrassenAgainstStandard(512, false);
		

	/* Time the standard Matrix Multiply vs N*/
		timeStandardMultiply();
	/*Time the Strassen Multiply vs N*/
		timeStrassenMultiply();	
	/* Blas Routine */
		timeBlas();
		
}



void testMatrixMultiply(){
	for(int p=3; p<=11; p++){
		int N = 1<<p;
		double* A = genRandomMat(N);
		double* B = genRandomMat(N);
		double start, time1, time2;
		
		start = omp_get_wtime();
		double* stdMatMul1 = mat_mul(N,A,B);
		time1 = omp_get_wtime()-start;
				
		start = omp_get_wtime();
		double* stdMatMul2 = mat_mul_menManage(N,A,B);
		time2 = omp_get_wtime()-start;
		
		printf("%d %g %g \n", N, time1, time2);
		
		delete A;
		delete B;
		delete stdMatMul2;
		delete stdMatMul1;
	}
}

void checkBLASStrassenAgainstStandard(int N, bool Print){
	double* A = genRandomMat(N);
	double* B = genRandomMat(N);
	double* C1 = mat_mul_menManage(N, A, B);
	double* C2 = writeStrassen2(N, A, B);
	double* C3 = BlasMultiply(N,A,B); 
	
	if(Print){
		printf("Matrix A \n");
		printMat(N,A);
		printf("Matrix B \n");
		printMat(N,B);
		printf("Standard Matrix Mulitply \n");
		printMat(N, C1);
		printf("Strassen Multiply Output \n");
		printMat(N,C2);
		printf("Blas dgemm Multiply Output \n");
		printMat(N,C3);
	}
	// check each element comparison
	int errS=0, errB=0;
	for(int i=0; i<N*N; i++){
		if(C1[i]!=C2[i]){
			errS+=1;
		}
		if(C1[i]!=C3[i]){
			errB+=1;
		}
	}
	
	if(errS==0) printf("N = %d; Standard and Strassen Match\n", N);
	else printf("N = %d; Standard and Strassen Does not Match!\n", N);
	
	if(errB==0) printf("N = %d; Standard and BLAS Match\n", N);
	else printf("N = %d; Standard and BLAS Does not Match!\n", N);
	
	delete A;
	delete B;
	delete C1;
	delete C2;
	delete C3;
}

void timeStandardMultiply(){
	for(int p=3; p<=12; p++){
		int N = 1<<p;
		double* A = genRandomMat(N);
		double* B = genRandomMat(N);
		double start, time1;
		
		start = omp_get_wtime();
		double* stdMatMul2 = mat_mul_menManage(N,A,B);
		time1 = omp_get_wtime()-start;
				
		printf("%d %g \n", N, time1);
		
		delete A;
		delete B;
		delete stdMatMul2;
	}
}

void timeStrassenMultiply(){
	for(int p=3; p<=12; p++){
		int N = 1<<p;
		double* A = genRandomMat(N);
		double* B = genRandomMat(N);
		double* C = new double[N*N];
		double time, start;
		
		start = omp_get_wtime();
		C = writeStrassen2(N, A, B);
		time = omp_get_wtime()-start;
		printf("%d %g \n", N, time);
		
		delete A;
		delete B;
		delete C;
	}
}

void timeBlas(){
	int rmax=1;

	for(int p=3; p<=12; p++){
		double time=0.0;
		int N = 1<<p;
		for(int r=0; r<rmax; r++){
			double* A = genRandomMat(N);
			double* B = genRandomMat(N);
			double* C = new double[N*N];
			double start=0;
			
			start = omp_get_wtime();
			C = BlasMultiply(N,A,B); 
			time += omp_get_wtime()-start;
			delete A;
			delete B;
			delete C;
		}
		printf("%d %g \n", N, time/(double)rmax);
	}
}



















