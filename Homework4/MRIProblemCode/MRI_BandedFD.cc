// Dean Hazineh 

#include "MRI_BandedFD.hh"
#include <cstdio>
#include <cmath>
#include <string>
#include <cstring>
#include <iostream>
#include <vector> //
#include <algorithm> // std::fill_n
#include <complex>   // std::complex
#include "file_output.hh"
#include <stdexcept>


// Define a wrapper for the banded lapack solver to be used by the class
typedef std::complex<double> dcomplex;
extern "C" { int zgbsv_(const int *n, const int *kl, const int *ku, const int *nrhs, dcomplex *ab, const int *ldab, int *ipiv, dcomplex *b, const int *ldb, int *info); }


MRI_BandedFD::MRI_BandedFD(int n_,double freq) : n(n_), nn(n_*n_), kl(n_), ku(n_), ldab(2*kl+ku+1), TB(new dcomplex[ldab*nn]), f(new dcomplex[nn]), ipiv(new int[nn]), diags(ldab*nn, 0.), h(1./(n_+1)), hsquared(h*h), ksquared(new dcomplex[nn]){
	start = diags.begin()+nn*kl;
	omegaSquared = pow(2.*M_PI*freq,2.); 
	fillKVectorVacuum(); // initializes k vector to vacuum conditions e_r = 1+0j
}


MRI_BandedFD::~MRI_BandedFD(){
	delete[] TB;
	delete[] f;
	delete[] ipiv;
	delete[] ksquared;
}


void MRI_BandedFD::fillKVectorVacuum(){
	for(int i=0; i<nn; i++) ksquared[i] = dcomplex(omegaSquared*mu_naught*epsilon_naught, 0.0);
}

void MRI_BandedFD::fillKVectorFromFile256(){
	if(n!=257) throw std::invalid_argument( "256 MRI Load requires n=257" );
	
	dcomplex* e_r = new dcomplex[257*257];
   FILE* file = fopen("mri_data_256.dat", "rb");
   fread(e_r, sizeof(dcomplex), 257*257, file);
   fclose(file);
    
   //for(int i=0; i<nn; i++) ksquared[i] = dcomplex(omegaSquared*mu_naught*epsilon_naught*e_r[i].real(), omegaSquared*mu_naught*epsilon_naught*e_r[i].imag());
   for(int i=0; i<nn; i++) ksquared[i] = dcomplex(omegaSquared*mu_naught*epsilon_naught*e_r[i].real(), omegaSquared*mu_naught*epsilon_naught*e_r[i].imag());
   saveComplexData(ksquared, "Out/KsquaredReal", "Out/KsquaredImag", "Out/KsquaredMag");
   delete[] e_r;	
}


void MRI_BandedFD::saveComplexData(dcomplex* complexdata, char* nameReal, char* nameImag, char* nameMag){	
	double realdata[nn];
	double imagdata[nn];
	double magdata[nn];
	for (int i=1; i<nn; i++){
		double thisreal = complexdata[i].real();
		double thisimag = complexdata[i].imag();
		realdata[i] = thisreal;
		imagdata[i] = thisimag;
		magdata[i] = sqrt(thisreal*thisreal + thisimag*thisimag);
	} 
	gnuplot_output(nameReal, realdata,n,n,ax,bx,ay,by); 
	gnuplot_output(nameImag, imagdata,n,n,ax,bx,ay,by); 
	gnuplot_output(nameMag, magdata,n,n,ax,bx,ay,by); 
}


void MRI_BandedFD::bandedSolveLapack(bool print){
	// Run LU solver with Lapack
   zgbsv_(&nn, &kl, &ku, &nrhs, TB, &ldab, ipiv, f, &nn, &info);

   // Process errors and Print if runs
	if (info < 0) {
		printf("\nInvalid argument at position %d\n", -info);
		}else if (info > 0) {
			printf("\nMatrix is singular at position %d\n", info);
      }else{
			if(print) saveComplexData(f, "Out/BandedSol_real", "Out/BandedSol_imag", "Out/BandedSol_mag");				
		}
}


void MRI_BandedFD::fillSourceTerm(bool print, int switchcase){
	for(int i=0; i<nn; i++) f[i] = dcomplex(0.0, 0.0);
	
	switch(switchcase){
	// Case 0 provides a Plus and Minus Step Function
		case 0:{
			int cc = floor(n/2);
			int sql = floor(n/3);
			for(int r=-sql; r<sql; r++){
				for(int c=-sql; c<sql; c++){
					if(r<0.) f[(cc+r)*n+(cc+c)] = dcomplex(1.0, 0.25);
					else f[(cc+r)*n+(cc+c)] = dcomplex(-1.0, -0.25);
				}	
			}
			break;
		}
		
	// Case 1 provides the dirac impulse at the center
		case 1:{
			int cc = floor(n/2);
			f[(cc)*n + cc] = dcomplex(1.0,0.0);
			break;
		}
	
	// Case 2 -- utilize the manufactured solution scenario from lectures
		case 2:{
			for(int r=0; r<n; r++) {
        		double y=(r+1)*h;
        		for(int c=0; c<n; c++) {
            	double x=(c+1)*h;
            	f[c+n*r]= dcomplex(-exp(x)*x*(-2-3*y+3*y*y+x*(2-y+y*y)), 0.0);
     			}
 			}
			break;
		}
	// Case 3 -- Place delta(x-(0.6,0.7)) source as required in MRI Problem
		case 3:{
			int row = floor(n*0.6);
			int col = floor(n*0.7);
			f[col+n*row] = dcomplex(1.0, 0.0);
			break;
		}
		
	}
	
	if(print) saveComplexData(f, "Out/BandedSource_real", "Out/BandedSource_imag", "Out/BandedSource_mag");
}


void MRI_BandedFD::fillTransferMatrix(bool print){
	 // Use diags and start iterator vector to aid in filling the banded matrix form for lapack

	 // start with main diagonal 
    //std::fill_n(start+nn*ku, nn, dcomplex(4./hsquared,0.)+ksquared); 
    int diagStart = nn*kl+nn*ku;
	 for(int i = 0; i<nn; i++) diags[diagStart+i] = dcomplex(4./hsquared + ksquared[i].real(), 0.+ksquared[i].imag());
    

    // Fill upper bands starting with the first upper through ku upper bandwidth
    for(int bandN = 1; bandN<=ku; bandN++){
    	int setVal = 0.0;
    	if(bandN==1||bandN == ku) setVal = -1./hsquared;   	
    	std::fill_n(start+nn*(ku-bandN)+bandN, nn-bandN, dcomplex(setVal, 0.0) );
    }
    // impose dirichlet boundary conditions right wall
    for(int i = n; i<nn; i+=n) *(start+nn*(ku-1) + i)=dcomplex(0.0,0.0);


    // Fill lower bands starting with first lower through kl
    for(int bandN = 1; bandN<=kl; bandN++){
    	int setVal = 0.0;
    	if(bandN==1 || bandN == kl) setVal = -1./hsquared;   	
    	std::fill_n(start + nn*(ku+bandN), nn-bandN, dcomplex(setVal, 0.0) );
    }
    // impose dirichlet boundary conditions left wall
    for(int i = n; i<nn; i+=n) *(start+nn*(ku+1) + (i-1))=dcomplex(0.0,0.0);
    
    
    // Transpose resulting diags holder to column-major order for LAPACK although we 
    // mentally wrote the stencil in row-major indexing!
    for (int r=0; r<ldab; r++)
        for (int c=0; c<nn; c++)
            TB[c*ldab+r] = diags[c+r*nn];
    
    
    // Print for checks and demonstration
    if(print){
		 // print the row major banded matrix
		 printf("Display the Row-Major Banded Matrix Form \n");
		 for(int r = 0; r<ldab; r++){
	 	 printf(r==0?"A=[":"  [");		
		 	for(int c = 0; c<nn; c++){
		 		dcomplex x = diags[c + r*nn];
	 			if(x.real()==0.0 && x.imag()==0.0)printf("  *   ");
		      else printf("%2.g%+2.gi ", x.real()*hsquared, x.imag()*hsquared);
		 	};
		 	puts("]");
		 };
		 puts("\n");
		 
		 // Print the matrix
	 	 printf("Display the FD Stencil Matrix used: \n");
		 for (int j=0; j<nn; j++) {
		     printf(j==0?"A=[":"  [");
		     for (int i=0; i<nn; i++) {
		         int k=i-j;
		         dcomplex x = k<-kl||k>ku?0:TB[kl+ku+i*ldab-k];
   	 			if(x.real()==0.0 && x.imag()==0.0)printf("  *   ");
		         else printf("%2.g%+2.gi ", x.real()*hsquared, x.imag()*hsquared);
		     }
		     puts("]");
		 }
	}	
}


void MRI_BandedFD::viewMMS_Example(){
	dcomplex* trueSol = new dcomplex[nn];
	dcomplex* error = new dcomplex[nn];
	for(int r=0; r<n; r++) {
		double y=(r+1)*h;
		for(int c=0; c<n; c++) {
			double x=(c+1)*h;
			double thisval = exp(x)*x*(1-x)*y*(1-y);
			trueSol[c+n*r] = dcomplex(thisval,0.0);
         error[c+n*r] = dcomplex(f[c+n*r].real()-thisval,f[c+n*r].imag());
        }
    }
	 saveComplexData(trueSol, "Out/mmsSolReal", "Out/mmsSolImag", "Out/mmsSolMag");
 	 saveComplexData(error, "Out/BandedMMSError_real", "tempOut", "tempOut");
    delete[] trueSol;
    delete[] error;
}











