// Dean Hazineh
// A really redundant and poorly coded FFTW (plus iterative) solver 
// as a minimum working implementation 
//--> I know its a mess but v1 was much cleaner, albeit wrong!

#ifndef MRI_DST_FFTW_V2_HH
#define MRI_DST_FFTW_V2_HH

#include "file_output.hh"
#include <fftw3.h>
#include <vector> 
#include <complex>// std::complex
#include <cmath>
typedef std::complex<double> dcomplex;

class MRI_DST_FFTW_v2{
	public:
		const int n; // length along 1D of grid not including zero boundary
		const int nn; // total length of vectorized grid space
		const double freq; // Used in defining wavevector K
		const double omegaSquared;
		const double h; // grid step size
		
		
		double* fre; // The discretized source term
		double* fco; // The discretized source term
		
		double* vre; // The discretized solution 
		double* vco; // The discretized solution 
		
		double* wre; // Holds the Fourier transform of the source
		double* wco; // Holds the Fourier transform of the source
			
		dcomplex knaughtsquared; // Contains the constant vacuum wavevector -- complex to match implementation
		dcomplex* ksquared; // Contains the spatially varying, complex vector
			
		MRI_DST_FFTW_v2(int n_, double freq);
		~MRI_DST_FFTW_v2();
		void fillSourceTerm(bool print, int switchcase);

		void solve(bool print);
		void saveDoubleData(double* savedata, char* nameReal);
		void viewMMS_Example();
		void IterSolve(bool print, int maxIter);
		void iterUpdateSource();
	   double computeIterativeSolverError();
	private:
		void fillKVectorVacuum();
		void fillKVectorFromFile256();
		double* const lam;
		fftw_plan plan_fwdRe;
		fftw_plan plan_fwdCo;
		fftw_plan plan_bckRe;
		fftw_plan plan_bckCo;

		const double ax=0,bx=1,ay=0,by=1;
		const double dx=(bx-ax)/(n+1),dy=(by-ay)/(n+1);
		const double mu_naught = 4.*M_PI*1e-7;
		const double epsilon_naught = 8.8542e-12;
};

#endif
