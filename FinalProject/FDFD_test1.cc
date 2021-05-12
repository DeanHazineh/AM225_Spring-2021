// Dean Hazineh
// AP225 Specific

#include "FDFD.hh"
#include <cstdio>

#include <cmath>
#include <math.h>
#include <complex>// std::complex
typedef std::complex<double> dcomplex;


int main(){

	int nx_=20, ny_=20, nz_=2, nPML_= 1;
	double dx_=1e-6, dy_=1e-6, dz_=1e-6, dPML_= 1e-6, omega_=.5e-6;
	FDFD_grid* MESolver_FDFD3D = new FDFD_grid(omega_, dx_, dy_, dz_, dPML_, nx_, ny_, nz_, nPML_);

	MESolver_FDFD3D->fillSpaceMaterial();
	MESolver_FDFD3D->viewEps();
	MESolver_FDFD3D->defineSourceTerm();
	MESolver_FDFD3D->viewJ();
	
	
	//MESolver_FDFD3D->AssembleAMatrixDoubleDemo();
	MESolver_FDFD3D->AssembleAMatrix();
	MESolver_FDFD3D->BICGSolve();
	MESolver_FDFD3D->	viewESol();
	
}
