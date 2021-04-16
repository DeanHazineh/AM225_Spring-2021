#ifndef MRI_BANDEDFD_HH
#define MRI_BANDEDFD_HH

#include <vector> 
#include <complex>// std::complex
#include <cmath>
typedef std::complex<double> dcomplex;

class MRI_BandedFD{
	public: 
		const int n; // length along 1D of grid not including zero boundary
		const int nn; // total length of vectorized grid space
		const int kl; // lower bandwidth
		const int ku; // upper bandwidth 
		const int ldab;
		const double h;
		const double hsquared;
		double omegaSquared;
		
		dcomplex* TB; // Common notation to have end B to denote banded matrix structure
		dcomplex* f; // Initialize source term and hold complex solution 
		dcomplex* ksquared; // contain the possibly complex wave-vector 

		MRI_BandedFD(int n_, double freq);
		~MRI_BandedFD();
		void fillKVectorVacuum();
		void fillKVectorFromFile256();
		void saveComplexData(dcomplex* complexdata, char* nameReal,char* nameImag, char* nameMag);
		void fillTransferMatrix(bool print);
		void fillSourceTerm(bool print, int switchcase);
		void bandedSolveLapack(bool print);
		void viewMMS_Example();

	private:
		// Initialize the internal parameters required for LAPACK Band Solver Below: 
		// Number of lower and upper diagonals
		// Number of source terms
		const int nrhs = 1;
		// Number of rows in banded storage
		// Pivot storage
		int* ipiv;
		// Error flag
		int info;	
		std::vector<dcomplex> diags;
		std::vector<dcomplex>::iterator start;
		
		// Initialize Parameters for outputting to file_output
		const double ax=0,bx=1,ay=0,by=1;
	   const double dx=(bx-ax)/(n+1),dy=(by-ay)/(n+1);
		const double mu_naught = 4.*M_PI*1e-7;
		const double epsilon_naught = 8.8542e-12;
};

#endif
