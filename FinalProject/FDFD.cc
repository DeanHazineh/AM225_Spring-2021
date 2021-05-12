// Dean Hazineh
// AP225 Specific
#include "FDFD.hh"
#include "omp.h"
#include <cstdio>
#include "file_output.hh"
#include <cmath>
#include <math.h>
#include <complex>// std::complex
typedef std::complex<double> dcomplex;
#include <stdexcept>

// *****************************************************************************
// ****************** Class Constructor ****************************************
//******************************************************************************
FDFD_grid::FDFD_grid(double omega_, double dx_, double dy_, double dz_, double dPML_, int nx_, int ny_, int nz_, int nPML_):
	omega(omega_),
	dx(dx_), dy(dy_), dz(dz_), dPML(dPML_),
	nx(nx_), ny(ny_), nz(nz_), nPML(nPML_),
	Lx(nx*dx), Ly(ny*dy), Lz(nz*dz), LPML(nPML*dPML),
	gridLen(3*(2*nPML+nx)*(2*nPML+ny)*(2*nPML+nz)), fieldOriginIdx(3*(2*nPML+nx)*(2*nPML+ny)*nPML+3*(2*nPML+nx)*nPML+3*nPML), 	
	zsliceGridSize((2*nPML+nx)*(2*nPML+ny)),
	e(new dcomplex[gridLen]), j(new dcomplex[gridLen]), 
	diag_mu(new dcomplex[gridLen]), diag_eps(new dcomplex[gridLen]),
	Ch_real(new double[gridLen*gridLen]), Ch_imag(new double[gridLen*gridLen]), Ce_real(new double[gridLen*gridLen]), Ce_imag(new double[gridLen*gridLen]),
	Areal(new double[gridLen*gridLen]), Aimag(new double[gridLen*gridLen]), m(4.),interm(new double[gridLen*gridLen]), interm2(new double[gridLen*gridLen]),
	A(new dcomplex[gridLen*gridLen])
	{
		sigma=(m+1.)*1.*log(0.75)/(2.*sqrt(munaught/epsilonnaught)*LPML*nPML);
		printf("Running FDFD_Test1\n");	
		printf("nx ny, nz: %d, %d, %d \n", nx, ny, nz );
		printf("Ni Nj, Nk: %d, %d, %d \n", (2*nPML+nx),(2*nPML+ny),(2*nPML+nz) );
		printf("Simulation Space Length Lx, Ly, Lz LPML: %g, %g, %g %g \n", Lx, Ly, Lz, LPML);		
		printf("Simulation grid discretizations dx, dy, dz dPML: %g, %g, %g %g \n", dx, dy, dz, dPML);		
		printf("Total Grid Values (Spatial x 3 Vector Components): %d \n", gridLen);
		printf("nth index of first nonPML origin point: %d \n", fieldOriginIdx);		
		printf("matrix A require size: %d \n", gridLen*gridLen);		
		printf("SigmaMax Setting: %f \n", sigma);
		
	}
	
FDFD_grid::~FDFD_grid(){
	delete[] A;
	delete[] interm;
	delete[] interm2;
	delete[] Aimag;
	delete[] Areal;
	delete[] Ce_imag;
	delete[] Ch_imag;
	delete[] Ce_real;
	delete[] Ch_real;
	delete[] diag_eps;
	delete[] diag_mu;
	delete[] j;
	delete[] e;
	}

// *****************************************************************************
// ****************** Initialize Problem Matrix Ax=B ***************************
//******************************************************************************
void FDFD_grid::AssembleAMatrixDoubleDemo(){
	// First create Ce (Electric Curl Operator Matrix)
	// NOTE: Ce here already includes multiplication by matrix Du^-1
	// create a pointer to each gridlen row
	double* cptr = Ce_real;
	dcomplex* useMu = diag_mu;
	for(int k=0; k<(nz+2*nPML); k++){
		for(int j=0; j<(ny+2*nPML); j++){
			for(int i=0; i<(nx+2*nPML); i++){
				
				//set the Curl for (x,i,j,k)
				if(j<(ny+2*nPML)-1) *(cptr+(2*nPML+nx)*3+2)    = 1./dy/useMu[0].real();
				 *(cptr+2)                                     = -1./dy/useMu[0].real();
				if(k<(nz+2*nPML)-1) *(cptr+zsliceGridSize*3+1) = -1./dz/useMu[0].real();
				*(cptr+1)                                      = 1./dz/useMu[0].real();
				
				cptr = cptr + 1*(gridLen)+1;
				// set the Curl for (y,i,j,k)
				if(k<(nz+2*nPML)-1) *(cptr+zsliceGridSize*3-1) = 1./dz/useMu[1].real();
				*(cptr-1)                                      = -1./dz/useMu[1].real();
				if(i<(nx+2*nPML)-1) *(cptr+1+3)                = -1./dx/useMu[1].real();
				*(cptr+1)                                      = 1./dx/useMu[1].real();

				cptr = cptr + 1*(gridLen)+1;			
				// set the curl for (z,i,j,k)	
				if(i<(nx+2*nPML)-1) *(cptr+3-1)                = 1./dx/useMu[2].real();
				*(cptr-1)                                      = -1./dx/useMu[2].real();
				if(j<(ny+2*nPML)-1) *(cptr+(2*nPML + nx)*3 -2) = -1./dy/useMu[2].real();
				*(cptr-2)                                      = 1./dy/useMu[2].real();	

		  		cptr = cptr + 1*(gridLen)+1;			
		  		useMu=useMu+3;
			}
		}
	}
	
	// Next create Ch (Magnetic Curl Operator Matrix) 
	// create a pointer to each gridlen row
	cptr = Ch_real;
	for(int k=0; k<(nz+2*nPML); k++){
		for(int j=0; j<(ny+2*nPML); j++){
			for(int i=0; i<(nx+2*nPML); i++){
			
				//set the Curl for (x,i,j,k)
				*(cptr+2)                          = 1./dy;
				if(j>0) *(cptr+2-(2*nPML+nx)*3)    = -1./dy;
				*(cptr+1)                          = -1./dz;
				if(k>0) *(cptr+1-zsliceGridSize*3) = 1./dz;
				
				cptr = cptr + 1*(gridLen)+1;
				// set the Curl for (y,i,j,k)
				*(cptr-1)                          = 1./dz;
				if(k>0) *(cptr-1-zsliceGridSize*3) = -1./dz;
				*(cptr+1)                          = -1./dx;
				if(i>0) *(cptr+1-3)                = 1./dx;
				
				cptr = cptr + 1*(gridLen)+1;			
				// set the curl for (z,i,j,k)	
				*(cptr-1)                          = 1./dx;
				if(i>0) *(cptr-1-3)                = -1./dx;
				*(cptr-2)                          = -1./dy;
				if(j>0) *(cptr-2-(nx+2*nPML)*3)    = 1./dy; 
				
		  		cptr = cptr + 1*(gridLen)+1;			
			}
		}
	}
	
	// Print Matrices to view in gnuplot as a check during dev.
	displayTransferMatrix(Ce_real, "Ce_real");
	displayTransferMatrix(Ch_real, "Ch_real");
	// Compute Ch*(Du^-1Ce) - omega^2De
	// output is initialized to omega^2De and will hold the output
	double* outptr = Areal;
	dcomplex* epsptr = diag_eps;
	for(int n=0; n<gridLen; n++){
		*(outptr) = omega*omega*epsptr[0].real();
		outptr = outptr+1+gridLen;
		epsptr=epsptr+1;
	}
	displayTransferMatrix(Areal, "DeReal");

	BlasMultiply(gridLen, Ch_real, Ce_real, Areal, 1.0, -1.0);
	displayTransferMatrix(Areal, "Areal");
}


void FDFD_grid::AssembleAMatrix(){
	// The matrcices will all in general be complex meaning we need to be careful 
	// The curl matrix will pick up complex values at the PML Boundaries!!!
	
	// First create Ce (Electric Curl Operator Matrix)
	// NOTE: Ce here already includes multiplication by matrix Du^-1
	double* cptr = Ce_real;
	double* cptri = Ce_imag;
	dcomplex* useMu = diag_mu;
	double sreal = 0., simag = 0., L = 0., sc = 0.;
	for(int k=0; k<(nz+2*nPML); k++){
		for(int j=0; j<(ny+2*nPML); j++){
			for(int i=0; i<(nx+2*nPML); i++){
				
				// Set Ce //
				//set the Curl for (x,i,j,k)
					if(i<nPML+nx && i>nPML){// if not in PML
						sreal = 1.0;
						simag = 0.0;					
					}else{
						if(i>nPML+nx) L = (i-nPML+nx)*LPML;
						else L=i*LPML;
						sc = sigma*pow(L/(nPML*LPML),m)/omega/epsilonnaught;
						sreal = 1./(1.+sc*sc);
						simag = sc/(1.+sc*sc);
					}
					//set real
					if(j<(ny+2*nPML)-1) *(cptr+(2*nPML+nx)*3+2)    = 1./dy/useMu[0].real()*sreal;
					*(cptr+2)                                     = -1./dy/useMu[0].real()*sreal;
					if(k<(nz+2*nPML)-1) *(cptr+zsliceGridSize*3+1) = -1./dz/useMu[0].real()*sreal;
					*(cptr+1)                                      = 1./dz/useMu[0].real()*sreal;
					//set imag
					if(j<(ny+2*nPML)-1) *(cptri+(2*nPML+nx)*3+2)    = 1./dy/useMu[0].imag()*simag;
					*(cptri+2)                                     = -1./dy/useMu[0].imag()*simag;
					if(k<(nz+2*nPML)-1) *(cptri+zsliceGridSize*3+1) = -1./dz/useMu[0].imag()*simag;
					*(cptri+1)                                      = 1./dz/useMu[0].imag()*simag;
				
						
				// set the Curl for (y,i,j,k)					
					cptr = cptr + 1*(gridLen)+1;
					cptri = cptri + 1*(gridLen)+1;
					if(j<nPML+ny && j>nPML){// if not in PML
						sreal = 1.0;
						simag = 0.0;					
					}else{
						if(j>nPML+ny) L = (j-nPML+ny)*LPML;
						else L=j*LPML;
						sc = sigma*pow(L/(nPML*LPML),m)/omega/epsilonnaught;
						sreal = 1./(1.+sc*sc);
						simag = sc/(1.+sc*sc);
					}
					// set real
					if(k<(nz+2*nPML)-1) *(cptr+zsliceGridSize*3-1) = 1./dz/useMu[1].real()*sreal;
					*(cptr-1)                                      = -1./dz/useMu[1].real()*sreal;
					if(i<(nx+2*nPML)-1) *(cptr+1+3)                = -1./dx/useMu[1].real()*sreal;
					*(cptr+1)                                      = 1./dx/useMu[1].real()*sreal;
					// set imag
					if(k<(nz+2*nPML)-1) *(cptri+zsliceGridSize*3-1) = 1./dz/useMu[1].imag()*simag;
					*(cptri-1)                                      = -1./dz/useMu[1].imag()*simag;
					if(i<(nx+2*nPML)-1) *(cptri+1+3)                = -1./dx/useMu[1].imag()*simag;
					*(cptri+1)                                      = 1./dx/useMu[1].imag()*simag;


				// set the curl for (z,i,j,k)	
					cptr = cptr + 1*(gridLen)+1;		
					cptri = cptri + 1*(gridLen)+1;
					if(k<nPML+nx && k>nPML){// if not in PML
						sreal = 1.0;
						simag = 0.0;					
					}else{
						if(k>nPML+nz) L = (k-nPML+nz)*LPML;
						else L=k*LPML;
						sc = sigma*pow(L/(nPML*LPML),m)/omega/epsilonnaught;
						sreal = 1./(1.+sc*sc);
						simag = sc/(1.+sc*sc);
					}
					// set real
					if(i<(nx+2*nPML)-1) *(cptr+3-1)                = 1./dx/useMu[2].real()*sreal;
					*(cptr-1)                                      = -1./dx/useMu[2].real()*sreal;
					if(j<(ny+2*nPML)-1) *(cptr+(2*nPML + nx)*3 -2) = -1./dy/useMu[2].real()*sreal;
					*(cptr-2)                                      = 1./dy/useMu[2].real()*sreal;	
					// set imag
					if(i<(nx+2*nPML)-1) *(cptri+3-1)                = 1./dx/useMu[2].imag()*simag;
					*(cptri-1)                                      = -1./dx/useMu[2].imag()*simag;
					if(j<(ny+2*nPML)-1) *(cptri+(2*nPML + nx)*3 -2) = -1./dy/useMu[2].imag()*simag;
					*(cptri-2)                                      = 1./dy/useMu[2].imag()*simag;	

					cptr = cptr + 1*(gridLen)+1;		
					cptri = cptri + 1*(gridLen)+1;
			  		useMu=useMu+3;
			}
		}
	}
	
	
	// Next create Ch (Magnetic Curl Operator Matrix) 
	// create a pointer to each gridlen row
	cptr = Ch_real;
	cptri = Ch_imag;
	sreal = 0., simag = 0., L = 0., sc = 0.;
	for(int k=0; k<(nz+2*nPML); k++){
		for(int j=0; j<(ny+2*nPML); j++){
			for(int i=0; i<(nx+2*nPML); i++){
				
				// Set Ch //
				//set the Curl for (x,i,j,k)
					if(i<nPML+nx && i>nPML){// if not in PML
						sreal = 1.0;
						simag = 0.0;					
					}else{
						if(i>nPML+nx) L = (i-nPML+nx)*LPML;
						else L=i*LPML;
						sc = sigma*pow(L/(nPML*LPML),m)/omega/epsilonnaught;
						sreal = 1./(1.+sc*sc);
						simag = sc/(1.+sc*sc);
					}
					//set real
					*(cptr+2)                          = 1./dy*sreal;
					if(j>0) *(cptr+2-(2*nPML+nx)*3)    = -1./dy*sreal;
					*(cptr+1)                          = -1./dz*sreal;
					if(k>0) *(cptr+1-zsliceGridSize*3) = 1./dz*sreal;
					//set imag
					*(cptri+2)                          = 1./dy*simag;
					if(j>0) *(cptri+2-(2*nPML+nx)*3)    = -1./dy*simag;
					*(cptri+1)                          = -1./dz*simag;
					if(k>0) *(cptri+1-zsliceGridSize*3) = 1./dz*simag;
				
				// set the Curl for (y,i,j,k)					
					cptr = cptr + 1*(gridLen)+1;
					cptri = cptri + 1*(gridLen)+1;
					if(j<nPML+nx && j>nPML){// if not in PML
						sreal = 1.0;
						simag = 0.0;					
					}else{
						if(j>nPML+ny) L = (j-nPML+ny)*LPML;
						else L=j*LPML;
						sc = sigma*pow(L/(nPML*LPML),m)/omega/epsilonnaught;
						sreal = 1./(1.+sc*sc);
						simag = sc/(1.+sc*sc);
					}
					// set real
					*(cptr-1)                          = 1./dz*sreal;
					if(k>0) *(cptr-1-zsliceGridSize*3) = -1./dz*sreal;
					*(cptr+1)                          = -1./dx*sreal;
					if(i>0) *(cptr+1-3)                = 1./dx*sreal;
					// set real
					*(cptri-1)                          = 1./dz*simag;
					if(k>0) *(cptri-1-zsliceGridSize*3) = -1./dz*simag;
					*(cptri+1)                          = -1./dx*simag;
					if(i>0) *(cptri+1-3)                = 1./dx*simag;
				

				// set the curl for (z,i,j,k)	
					cptr = cptr + 1*(gridLen)+1;		
					cptri = cptri + 1*(gridLen)+1;
					if(k<nPML+nx && k>nPML){// if not in PML
						sreal = 1.0;
						simag = 0.0;					
					}else{
						if(k>nPML+nz) L = (k-nPML+nz)*LPML;
						else L=k*LPML;
						sc = sigma*pow(L/(nPML*LPML),m)/omega/epsilonnaught;
						sreal = 1./(1.+sc*sc);
						simag = sc/(1.+sc*sc);
					}
					// set real
					*(cptr-1)                          = 1./dx*sreal;
					if(i>0) *(cptr-1-3)                = -1./dx*sreal;
					*(cptr-2)                          = -1./dy*sreal;
					if(j>0) *(cptr-2-(nx+2*nPML)*3)    = 1./dy*sreal; 
					// set imag
					*(cptri-1)                          = 1./dx*simag;
					if(i>0) *(cptri-1-3)                = -1./dx*simag;
					*(cptri-2)                          = -1./dy*simag;
					if(j>0) *(cptri-2-(nx+2*nPML)*3)    = 1./dy*simag; 
				
		  		cptr = cptr + 1*(gridLen)+1;		
				cptri = cptri + 1*(gridLen)+1;	
			}
		}
	}
	
	displayTransferMatrix(Ce_real, "Ce_real");
	displayTransferMatrix(Ch_real, "Ch_real");
	displayTransferMatrix(Ce_imag, "Ce_imag");
	displayTransferMatrix(Ch_imag, "Ch_imag");
	
	// Compute Ch*(Du^-1Ce) - omega^2De
	// A is initialized to omega^2De and will hold the output
	double* outptr = Areal;
	double* outptri = Aimag;
	dcomplex* epsptr = diag_eps;
	for(int n=0; n<gridLen; n++){
		*(outptr) = omega*omega*epsptr[0].real();
		*(outptri) = omega*omega*epsptr[0].imag();
		outptr = outptr+1+gridLen;
		outptri = outptri+1+gridLen;
		epsptr=epsptr+1;
	}
	displayTransferMatrix(Areal, "DeReal");
	displayTransferMatrix(Aimag, "Deimag");

	// Get the real part of the transfer Matrix
	BlasMultiply(gridLen, Ch_imag, Ce_imag, interm, 1.0, 0.0);	
	BlasMultiply(gridLen, Ch_real, Ce_real, interm, 1.0, -1.0);
	matrixAdd(gridLen*gridLen, interm, Areal, Areal, -1.0);
	
	// Get the imaginary part of the transfer Matrix
	BlasMultiply(gridLen, Ch_real, Ce_imag, interm2, 1.0, 0.0);	
	BlasMultiply(gridLen, Ch_imag, Ce_real, interm2, 1.0, 1.0);
	matrixAdd(gridLen*gridLen, interm2, Aimag, Aimag, -1.0);

	displayTransferMatrix(Areal, "Areal");
	displayTransferMatrix(Aimag, "Aimag");
	
	// From the dcomplex A matrix from the A real and A imaginary
	#pragma omp parallel for
	for(int n=0; n<gridLen*gridLen; n++){
		A[n] = dcomplex(Areal[n], Aimag[n]);
	}
}


// *****************************************************************************
// ****************** Initialize Source and Space Distribution Grid ************
//******************************************************************************

// NOTE: I know the formulation dictates to introduce source J which also includes definition of M
// from the equations written but I have no idea what these source terms look like or how to write them down
// for different scenarios. I just create something artificial source until I get a better idea
// Note, the factor -i\omega is included in J!
void FDFD_grid::defineSourceTerm(){
	dcomplex* jp = j;
	dcomplex zeroTerm = dcomplex(0.,0.);

	
	for(int k=0; k<(nz+2*nPML); k++){
		double zval = k*dz;
		for(int j=0; j<(ny+2*nPML); j++){
			double yval = j*dy;
			for(int i=0; i<(nx+2*nPML); i++){
				double xval = i*dx;	
				
				// Conditional Fill of Source Term
				//dcomplex sourceTerm = dcomplex(exp(-(xval*xval+yval*yval)/(2.*sphereRad*sphereRad)),exp(-(xval*xval+yval*yval)/(2.*sphereRad*sphereRad)));	
				dcomplex sourceTerm = dcomplex(10., 1./omega);
				if( zval>=LPML && zval < LPML+Lz){
					if(xval>=LPML && xval < LPML+Lx && yval>LPML && yval<=(LPML+Ly/7.)){
						*(jp)=sourceTerm;	
						*(jp+1)=sourceTerm;	
						*(jp+2)=sourceTerm;	
					}										
				}
				else{
					*(jp)=zeroTerm;
					*(jp+1)=zeroTerm;
					*(jp+2)=zeroTerm;
				}
				jp+=3;
			}
		}
	}
}

void FDFD_grid::fillSpaceMaterial(){
/*Note, mu and epsilon matrices should be converted to just double in the future and not complex */
	// Define electric permeability values
	double epsilon_relative = 2.;
	dcomplex freespace = dcomplex(epsilonnaught,epsilonnaught);
	dcomplex indexMaterial = dcomplex(epsilonnaught*epsilon_relative, epsilonnaught*epsilon_relative);
	// define Magnetic permebaility values	
	double mu_relative = 1.;
	dcomplex muInvariant = dcomplex(munaught,munaught);		
	
	// FIRST, initialize the entire epsilon matrix to free space 
	dcomplex* eps = diag_eps;
	dcomplex* mu = diag_mu;
	for(int k=0; k<(nz+2*nPML); k++){
		for(int j=0; j<(ny+2*nPML); j++){
			for(int i=0; i<(nx+2*nPML); i++){
				// set vacuum electric permittivity
				*(eps)= freespace;
				*(eps+1)= freespace;
				*(eps+2)= freespace;
				eps+=3;
				// set uniform magnetic permeability
				*(mu)=muInvariant;
				*(mu+1)=muInvariant;
				*(mu+2)=muInvariant;
				mu+=3;
			}
		}
	}
	
	// Add conditional material filler
	eps = diag_eps + (nPML)*3*zsliceGridSize;;
	double cx = (Lx+2*LPML)/2., cy=(Ly+2*LPML)/2., cz=(Lz+2*LPML)/2.;
	for(int k=nPML; k<nz+nPML; k++){
		double zval = k*dz - cz;
		for(int j=0; j<(ny+2*nPML); j++){
			double yval = j*dy - cy;
			for(int i=0; i<(nx+2*nPML); i++){
				double xval = i*dx - cx;
				
				if(xval*xval + yval*yval + zval*zval<=sphereRad*sphereRad){
					*(eps)= indexMaterial;
					*(eps+1)= indexMaterial;
					*(eps+2)= indexMaterial;
				}else{
					*(eps)= freespace;
					*(eps+1)= freespace;
					*(eps+2)= freespace;
				}
				eps+=3;
			}
		}
	}

}

// *****************************************************************************
// ****************** Try to write the BICG Solver   ***************************
//******************************************************************************
//https://en.wikipedia.org/wiki/Biconjugate_gradient_method
void FDFD_grid::BICGSolve(){
	// set up required constants
	dcomplex* minusOne = new dcomplex(-1.0, -1.0);
	dcomplex* one = new dcomplex(1.0, 1.0);	
	dcomplex* zero = new dcomplex(0.0, 0.0);
	int none = 1;
	
	dcomplex* r = new dcomplex[gridLen];
	dcomplex* rhat = new dcomplex[gridLen];
	dcomplex* rp = new dcomplex[gridLen];
	dcomplex* rphat = new dcomplex[gridLen];
	dcomplex* p = new dcomplex[gridLen];
	dcomplex* phat = new dcomplex[gridLen];
	dcomplex* x = new dcomplex[gridLen];
	dcomplex* xhat = new dcomplex[gridLen];
	dcomplex alpha = dcomplex(0.0,0.0);	
	dcomplex* nu = new dcomplex[gridLen];
	
	// set initial guess e to zeros
	for(int n=0; n<gridLen; n++){
		x[n] = dcomplex(0.0, 0.0);
	}
	ccopy_(&gridLen, x, &none, xhat, &none);	
	
	//RUN
	// implement ro <- b-Ax_o
		ccopy_(&gridLen, j, &none, r, &none);	
		// y = alpha*trans(A)*x+beta*y
		//     transA   nrowA     ncolA     alpha     A  LDA       x  incx   beta Y  incy
		cgemv_(&transn, &gridLen, &gridLen, minusOne, A, &gridLen, x, &none, one, r, &none);
	
	// implement rohat <- bhat - xhat*A 
		complexConjugateVector(gridLen, j, rhat);
		// C = alpha*op(A)*op(B) + beta*C
		//     Trans(A) trans(B) nrowA  ncolsB    ncolsA    Alpha     A  LDA    B  LDB       Beta C           LDC 
		cgemm_(&transt, &transn, &none, &gridLen, &gridLen, minusOne, xhat, &gridLen, A, &gridLen, one, rhat, &none);     
	
	// set po<-ro
 	// set pohat <- rhat
	 	ccopy_(&gridLen, r, &none, p, &none);
	 	ccopy_(&gridLen, rhat, &none, phat, &none);
 	
 	// RUN KRYLOV LOOP
 	dcomplex rsq = dcomplex(0.0,0.0);
 	dcomplex rhatsq = dcomplex(0.0,0.0);
 	for(int iter=0; iter<20; iter++){
	 rsq = cdotc_(&gridLen, r, &none, r, &none);
	 rhatsq = cdotc_(&gridLen, rhat, &none, rhat, &none);
    printf("iter %d, R: %g, Rhat: %g \n", iter, rsq.real(), rhatsq.real());

		// compute alpha <- rhat*r/(phat*A*p)
			cgemv_(&transn, &gridLen, &gridLen, one, A, &gridLen, p, &none, zero, nu, &none);
			alpha = cdotu_(&gridLen, rhat, &none, r, &none)/cdotu_(&gridLen, phat, &none, nu, &none);

			//dcomplex test = cdotu_(&gridLen, rhat, &none, r, &none);
			//printf("Test: %g %g\n", test.real(), test.imag());
			
		// xk+1 <- xk + alphak*pk
			//     lengthn   alpha,  X, incx, Y incy 
			caxpy_(&gridLen, &alpha, p, &none, x, &none );
		// xhatk+1 <- xhatk + alphak*phatk
			caxpy_(&gridLen, &alpha, phat, &none, xhat, &none );

		// rk+1 <- rk - alphak * Apk
			alpha = dcomplex(-1.*alpha.real(), -1.*alpha.imag());
		 	ccopy_(&gridLen, r, &none, rp, &none);
			cgemv_(&transn, &gridLen, &gridLen, &alpha, A, &gridLen, p, &none, one, rp, &none);
		// rhatk+1 <- rhat - alpha * p*A
 		 	ccopy_(&gridLen, rhat, &none, rphat, &none);
			cgemm_(&transt, &transn, &none, &gridLen, &gridLen, &alpha, phat, &gridLen, A, &gridLen, one, rphat, &none);
		
			// compute Beta_k <- (rphat dot rp/ rhat dot r)
			alpha = cdotu_(&gridLen, rphat, &none, rp, &none)/cdotu_(&gridLen, rhat, &none, r, &none);
			
			// pk+1 <- rp + beta*pk
			//     lengthn   alpha,  X, incx, Y incy 
		 	ccopy_(&gridLen, rp, &none, nu, &none);
			caxpy_(&gridLen, &alpha, p, &none, nu, &none );
		 	ccopy_(&gridLen, nu, &none, p, &none);
		 	ccopy_(&gridLen, rp, &none, r, &none);
		// phat <- rphat + beta*phat
		 	ccopy_(&gridLen, rphat, &none, nu, &none);
			caxpy_(&gridLen, &alpha, phat, &none, nu, &none );
		 	ccopy_(&gridLen, nu, &none, phat, &none);
		 	ccopy_(&gridLen, rphat, &none, rhat, &none);
 	}
 	if(rsq.real() < rhatsq.real()) ccopy_(&gridLen, x, &none, e, &none);
	else ccopy_(&gridLen, xhat, &none, e, &none);

 	delete[] nu; 	
 	delete[] xhat;
	delete[] x; 
	delete[] phat;
	delete[] p;
 	delete[] rphat;
 	delete[] rp;
	delete[] r;
	delete[] rhat;
}


// *****************************************************************************
// ********************** PRINT AND OUTPUT UTILITIES ***************************
//******************************************************************************
void FDFD_grid::viewEps(){
	output3D(diag_eps, "EpsilonMatrix");
}

void FDFD_grid::viewJ(){
	output3D(j, "JMatrix");
}

void FDFD_grid::viewESol(){
	output3D(e, "ESol");
}

void FDFD_grid::output3D(dcomplex* dataVector, const char* fileOutput){
	// To save space, just going to take the vector magnitude at (w,i,j,k) points
	// for each z slice and output it as a binary file excluding 
	#pragma omp parallel for
	for(int k=0; k<(nz+2*nPML); k++){
		double vectorMagnitude[zsliceGridSize];
	   char filename[50];
      sprintf(filename,"out/%s%01d",fileOutput,k);

		dcomplex* dp = dataVector + k*zsliceGridSize*3;
      for(int j=0; j<(ny+2*nPML); j++){
			for(int i=0; i<(nx+2*nPML); i++){
				vectorMagnitude[i+j*(nx+2*nPML)] =  sqrt( 
					(dp[0].real()*dp[0].real() + dp[0].imag()*dp[0].imag()) + 
					(dp[1].real()*dp[1].real() + dp[1].imag()*dp[1].imag()) + 
					(dp[2].real()*dp[2].real() + dp[2].imag()*dp[2].imag())
				);
			dp=dp+3;
			}
		}
		gnuplot_output(filename, vectorMagnitude, nx+2*nPML, ny+2*nPML, 0.,LPML*2+Lx, 0, LPML*2+Ly);
	}
	
}

void FDFD_grid::displayTransferMatrix(double* matrix, const char* fileOutput){
	char filename[50];
   sprintf(filename,"out/%s",fileOutput);
	gnuplot_output(filename, matrix, gridLen, gridLen, 0., gridLen, 0., gridLen);
}

void FDFD_grid::BlasMultiply(int N, double* A, double* B, double* C, double alpha, double beta){
  // Implement matrix multiplication with Blas
  // Due to C++/Fortran calling conventions, all entries are passed as pointers, and there is an
  // underscore after the name.
	char trans='n';
	dgemm_(&trans,&trans,&N,&N,&N,&alpha, A,&N,B,&N,&beta,C,&N);
}

//void FDFD_grid::BLAS_cgemv(const char trans, int m, int n, dcomplex* alpha, dcomplex *a, int lda, dcomplex* x, int incx,
//	dcomplex* beta, dcomplex *y, int incy){
//	// implement alpha*A*x + beta*y
//	cgemv_(&trans, &m, &n, alpha, a, &lda, x, &incx, beta, y, &incy);
//}

void FDFD_grid::matrixAdd(int N, double* A, double* B, double* C, double alpha){
	#pragma omp parallel for
	for(int n=0; n<N; n++){
		C[n] = A[n] + alpha*B[n];
	}
}

void FDFD_grid::complexConjugateVector(int N, dcomplex* a, dcomplex* b){
	#pragma omp parallel for
	for(int n =0; n<N; n++){
		b[n] = dcomplex(a[n].real(), -1.0*a[n].imag());
	}
}





























