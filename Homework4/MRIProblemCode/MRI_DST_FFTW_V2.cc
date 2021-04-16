// Made copy of the poisson_fft class from Rycroft 225 examples
// and inserted modifications and changes here

#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector> 
#include <complex>// std::complex
typedef std::complex<double> dcomplex;
#include "MRI_DST_FFTW_V2.hh"
#include "MRI_BandedFD.hh" // used to get reference solution for iterative error comparison

MRI_DST_FFTW_v2::MRI_DST_FFTW_v2(int n_, double freq) : 
	 freq(freq), n(n_), nn(n_*n_), h(1./(n+1.)), omegaSquared(pow(2.*M_PI*freq,2.)),
    ksquared(new dcomplex[nn]), lam(new double[n]), 
 	 fre(fftw_alloc_real(nn)), fco(fftw_alloc_real(nn)), 
 	 vre(fftw_alloc_real(nn)), vco(fftw_alloc_real(nn)), 
 	 wre(fftw_alloc_real(nn)), wco(fftw_alloc_real(nn)),
    plan_fwdRe(fftw_plan_r2r_2d(n,n,fre,wre,FFTW_RODFT00,FFTW_RODFT00,FFTW_MEASURE)),
    plan_fwdCo(fftw_plan_r2r_2d(n,n,fco,wco,FFTW_RODFT00,FFTW_RODFT00,FFTW_MEASURE)),
    plan_bckRe(fftw_plan_r2r_2d(n,n,wre,vre,FFTW_RODFT00,FFTW_RODFT00,FFTW_MEASURE)),
    plan_bckCo(fftw_plan_r2r_2d(n,n,wco,vco,FFTW_RODFT00,FFTW_RODFT00,FFTW_MEASURE))
    {
	 
	 // Define the constant k0^2 which is used in regular solver
	 knaughtsquared = dcomplex(omegaSquared*mu_naught*epsilon_naught, 0.0);
    
    // Initialize the table of eigenvalues for laplacian on sin terms
    //for(int j=0;j<n;j++) lam[j]=2*(1-cos(fac*(j+1)));
    const double fac = M_PI/(n+1);
    for(int j=0;j<n;j++) lam[j]=fac*fac*(j+1)*(j+1);
}

MRI_DST_FFTW_v2::~MRI_DST_FFTW_v2(){
    fftw_destroy_plan(plan_bckRe);
    fftw_destroy_plan(plan_bckCo);
    fftw_destroy_plan(plan_fwdRe);
    fftw_destroy_plan(plan_fwdCo);
    delete[] lam;
    delete[] ksquared;
    fftw_free(wco);
    fftw_free(wre);
    fftw_free(vco);
    fftw_free(vre);
    fftw_free(fre);
    fftw_free(fco);
}


void MRI_DST_FFTW_v2::fillKVectorFromFile256(){
/* Loads the MRI256 file which contains epsilon_relative across 256x256 grid */
	if(n!=257) throw std::invalid_argument( "256 MRI Load requires n=257" );
	
	dcomplex* e_r = new dcomplex[257*257];
   FILE* file = fopen("mri_data_256.dat", "rb");
   fread(e_r, sizeof(dcomplex), 257*257, file);
   fclose(file);
    
   // NOTE: To actual solve this with the iterative solver, it is useful to modify the matrix definition so 
   // ksquared = ksquared-knaughtsquared at each spatial location!
   for(int i=0; i<nn; i++) ksquared[i] = dcomplex(omegaSquared*mu_naught*epsilon_naught*(e_r[i].real()-1.0), omegaSquared*mu_naught*epsilon_naught*e_r[i].imag());
   delete[] e_r;	
}

void MRI_DST_FFTW_v2::saveDoubleData(double* savedata, char* nameReal){	
	gnuplot_output(nameReal,savedata,n,n,ax,bx,ay,by); 
}

void MRI_DST_FFTW_v2::fillSourceTerm(bool print, int switchcase){
/* This function must be called to intitialize the source term.
Boolean flag allows for printing the source term at the end */
	for(int i=0; i<nn; i++){
		fre[i] = 0.0;
		fco[i] = 0.0;
	}
	
	switch(switchcase){
		case 0:{
		// type 1 test // Not used anymore!
			break;}
		
		// Case 1 provides the dirac impulse at the center
		case 1:{
			int cc = floor(n/2);
			fre[(cc)*n + cc] = 1.0;
			fco[(cc)*n + cc] = 0.0;
			break;}	
			
		//Case 2 -- utilize the manufactured solution scenario from lectures
		case 2:{
			for(int r=0; r<n; r++){
        		double y=(r+1)*h;
        		for(int c=0; c<n; c++){
            	double x=(c+1)*h;
            	fre[c+n*r]= -exp(x)*x*(-2-3*y+3*y*y+x*(2-y+y*y));
            	fco[c+n*r]= 0.0;
     			}
 			}
			break;}
			
		// Case 3-- MRI Probe Problem
		case 3:{
			int row = floor(n*0.6);
			int col = floor(n*0.7);
			fre[col+n*row] = 1.0;
			fco[col+n*row] = 0.0;
			break;}
	}
	
	if(print){
		saveDoubleData(fre,"Out/FFTWSource_real");
		saveDoubleData(fco,"Out/FFTWSource_imag");
	}
}


/** Solves the linear system using the fast Fourier transform. 
This solve call only works for constant k!*/
void MRI_DST_FFTW_v2::solve(bool print){
	const double nor=1./(2*(n+1)),fac=h*h*nor*nor;
	fftw_execute(plan_fwdRe);
	fftw_execute(plan_fwdCo);

	for(int j=0;j<n;j++) for(int i=0;i<n;i++) {
	  wre[i+n*j]*=1./(lam[i]/fac + lam[j]/fac + knaughtsquared.real());
	  wco[i+n*j]*=1./(lam[i]/fac + lam[j]/fac + knaughtsquared.imag());
	}

	fftw_execute(plan_bckRe);
	fftw_execute(plan_bckCo);
	
	if(print){
	saveDoubleData(vre, "Out/FFTWSol_real");
	saveDoubleData(vco, "Out/FFTWSol_imag");
	}
}


void MRI_DST_FFTW_v2::IterSolve(bool print, int maxIter){
/* This function manages the iterative solver which allows for spatially
varying k vector */
   // fill spatially varying ksquared vector via loading MRI Data file
   fillKVectorFromFile256();
   
	for(int iter=0; iter<maxIter; iter++){
		// solve current iteration
		solve(false);
		
		// update the source terminal
		iterUpdateSource();
		
		// compute the error 
		printf("%d %g \n", iter, computeIterativeSolverError());
	}
	
	if(print){
    	saveDoubleData(vre, "Out/FFTWSol_real");
    	saveDoubleData(vco, "Out/FFTWSol_imag");
    }
}

void MRI_DST_FFTW_v2::iterUpdateSource(){
	// NOte: see definition of ksquared
	// Ksquaredvector = ksquared-knaughtsquared already subtracted
	
	fillSourceTerm(false, 3); // reset the source to regular v then add the modification to it 
	for(int i=0; i<nn; i++){
		fre[i] = fre[i] - (ksquared[i].real()*vre[i]-ksquared[i].imag()*vco[i]);
		fco[i] = fco[i] - (ksquared[i].real()*vco[i]+ksquared[i].imag()*vre[i]); 
	}		
}

double MRI_DST_FFTW_v2::computeIterativeSolverError(){
	MRI_BandedFD* LUBSol = new MRI_BandedFD(257, freq);
	LUBSol->fillKVectorFromFile256();
	LUBSol->fillTransferMatrix(false);
	LUBSol->fillSourceTerm(true, 3);
	LUBSol->bandedSolveLapack(true);
	
	double error1 = 0, error2 = 0, netError = 0;
	for(int i = 0; i<nn; i++){
		error1 = LUBSol->f[i].real()-vre[i];
		error2 = LUBSol->f[i].imag()-vco[i];
		netError+= (error1*error1 + error2*error2);
	}
	
	delete LUBSol;
	return netError; 
}


void MRI_DST_FFTW_v2::viewMMS_Example(){
	double* trueSol = new double[nn];
	double* error = new double[nn];
	for(int r=0; r<n; r++) {
		double y=(r+1)*h;
		for(int c=0; c<n; c++) {
			double x=(c+1)*h;
			double thisval = exp(x)*x*(1-x)*y*(1-y);
			trueSol[c+n*r] = thisval;
         error[c+n*r] = vre[c+n*r] - thisval;
        }
    }
	 saveDoubleData(trueSol,"Out/mmsSolReal");
 	 saveDoubleData(error,"Out/FFTWMMSError_real");
    delete[] trueSol;
    delete[] error;
}

















    
    

