/*
Dean Hazineh
Spring 2021
Code for Hw 3 Problem 2
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "rbfAnalysis.hh"

void PrintedDemoOfRBFSolver();
void computeMatrixCounts(int k, int mode, int hilbertM);
void viewRatioSorted();
void testHilbertArrangement();
void viewRatioHilbert();
void testHilbertArrangementSelect(int k, int m);
void rbfTimePCG();

int main(){
	/* Do a quick demo of the rbf code */
 	//PrintedDemoOfRBFSolver();

	/* print Hilbert sorted Data to plot in Gnuplot */
	//testHilbertArrangement();
	
	/* Generate Ratio Curves unsorted vs sorted */
	//viewRatioSorted();
	/* Generate Ratio using Hilbert Curver Sorting */
	//viewRatioHilbert();
	
	/* view selected Hilbert Curves */
	//testHilbertArrangementSelect(100000, 18);

	/* run time comparison of hilbert vs y sort */
	 rbfTimePCG();	
}



void PrintedDemoOfRBFSolver(){
	// run a full demo of the RBF Solvers with the settings
	// k=15, bls = 5, and using rycroft sorted initial coordinates
	 int k = 15;
	 rbf_test r(k,2);
	 r.initPoints(1,0);
	 
	 r.demo_LapackdenseSolver();
	 int Tcount = r.countDenseT();
	 printf("Number of Elements in Matrix: %d \n", Tcount);
	 printf("\n");
	 
	 r.demo_JacobiPCG(5);
	 int Pcount = r.countJacobiP(false, 5);
	 double ratio = (double) Pcount/ (double)Tcount;
  	 printf("k: %d Ratio P/T:%g Tcount: %d Pcount: %d \n", k, ratio, Tcount, Pcount);
	 printf("\n");	 
}

void computeMatrixCounts(int k, int mode, int hilbertM){
	rbf_test r(k,2);
	r.initPoints(mode, hilbertM);
   int Tcount = r.countDenseT();
   // use sqrt(k) as block size 
   int Pcount = r.countJacobiP(false, sqrt(k));
   double ratio = (double) Pcount/ (double)Tcount;
   printf("%d %g %d %d \n", k, ratio, Tcount, Pcount);
}

void viewRatioSorted(){
	for(int mode=0; mode<=1; mode++){//iterate unsorted and 1D sorted
		for(int k=10;k<10000;k+=k>>3){
			computeMatrixCounts(k,mode,0);
		}printf("\n");
	}
}

void testHilbertArrangement(){
	for(int M=0; M<=6; M++){
		rbf_test r(10000,2);
		r.init_Hilbert(M);
		r.print_HilbertData();
		printf("\n");
	}
}

void viewRatioHilbert(){
	int mode = 2;
	//for(int m=0; m<10; m+=1){
	for(int m=10; m<21; m+=1){
		for(int k=10;k<10000;k+=k>>2){
			computeMatrixCounts(k,mode,m);
		}printf("\n");
	}	
}

void testHilbertArrangementSelect(int k, int m){
	rbf_test r(k,2);
	r.init_Hilbert(m);
	r.print_HilbertData();
	printf("\n"); puts("");
}

void rbfTimePCG(){
	for(int mode=0; mode<=2; mode++){ // loop mode 0 = unsorted 1=sorted y and mode 2=sort Hilbert
		for(int k=10;k<10000;k+=k>>3){
			rbf_test r(k,2);
			r.initPoints(mode, 10);//use hilbert order 10 when hilbert mode
			double tcg = r.timeSolve_JacobiPCG(mode);
			printf("%d %g \n", k, tcg);
		}puts("\n");puts("");
	}
}




