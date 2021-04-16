#include "MRI_BandedFD.hh"
//#include "MRI_DST_FFTW.hh"
#include "MRI_DST_FFTW_V2.hh"
#include "omp.h" // used for timing function

void runExample_PrintFDStencil();
void runExample_MMSPoisson();
double runBandedGeneral(int gridn, double freq, bool save, int sourcetype);
double runDSTGeneral(int gridn, double freq, bool save, int sourcetype);
void runExample_delta();
void runTimedDeltaFunction();
void runMRI256(double freq);
void runIterativeFFT(int steps, double freq);


int main() {
	/* Print the FD stencil used in Banded Solver */
	//runExample_PrintFDStencil();
	
	/* Test against the MMS Source from Lecture*/
	//runExample_MMSPoisson();
	
	/* Run the dirac impulse test with non-zero k vector constant*/
	//runExample_delta();
	
	/* Time operations */
	//runTimedDeltaFunction();
	
	/* Run 256 MRI Example */
	//runMRI256(298.3e6);
	
	/* Run Lapack Iterative test*/
	runIterativeFFT(5, 298.3e6);

}


void runExample_PrintFDStencil(){
	MRI_BandedFD* test = new MRI_BandedFD(4, 0.0);
	test->fillTransferMatrix(true);
	delete test;
}

void runExample_MMSPoisson(){
	double freq = 0.0;
	int sourcetype = 2;
	int gridn = 250;
	
	MRI_BandedFD* test = new MRI_BandedFD(gridn, freq);
	test->fillTransferMatrix(false);
	test->fillSourceTerm(true,sourcetype);
	test->bandedSolveLapack(true);
	test->viewMMS_Example();
	delete test;
	
	MRI_DST_FFTW_v2* test2 = new MRI_DST_FFTW_v2(gridn, freq);
	test2->fillSourceTerm(true, sourcetype);
	test2->solve(true);
	test2->viewMMS_Example();
	delete test2;
}

void runExample_delta(){
	double freq = 21.3e6;
	int sourcetype = 1; // centered real delta function
	int gridn = 256;
	
	runBandedGeneral(gridn, freq, true, sourcetype);
	runDSTGeneral(gridn, freq, true, sourcetype);
}

double runBandedGeneral(int gridn, double freq, bool save, int sourcetype){
	double start;
	start = omp_get_wtime();
	MRI_BandedFD* test = new MRI_BandedFD(gridn, freq);
	test->fillTransferMatrix(false);
	test->fillSourceTerm(save, sourcetype);
	test->bandedSolveLapack(save);
	delete test;
	return omp_get_wtime()-start;
}

double runDSTGeneral(int gridn, double freq, bool save, int sourcetype){
	double start;
	start = omp_get_wtime();
	MRI_DST_FFTW_v2* test2 = new MRI_DST_FFTW_v2(gridn, freq);
	test2->fillSourceTerm(save, sourcetype);
	test2->solve(save);
	return omp_get_wtime()-start;
}


void runTimedDeltaFunction(){
	double freq = 21.3e6;
	bool save = false;
	int sourcetype = 1;
	// Time DST Solver
	for(int n=10;n<=1024;n+=n>>1) {
		double timeDST = runDSTGeneral(n, freq, save, sourcetype);
		printf("%d %g \n", n, timeDST);
	}
	puts("\n");
	// Time Banded FD Solver over smaller range
	for(int n=10;n<=550;n+=n>>1) {
		double timeBanded = runBandedGeneral(n, freq, save, sourcetype);
		printf("%d %g \n", n, timeBanded);
	}
}

void runMRI256(double freq){
	MRI_BandedFD* test = new MRI_BandedFD(257, freq);
	test->fillKVectorFromFile256();
	test->fillTransferMatrix(false);
	test->fillSourceTerm(true, 3);
	test->bandedSolveLapack(true);
	delete test;
}

void runIterativeFFT(int steps, double freq){
	MRI_DST_FFTW_v2* test = new MRI_DST_FFTW_v2(257, freq);
	test->fillSourceTerm(true, 3);
	test->IterSolve(true, steps);
	delete test;
}




