#ifndef KURAMOTOMODEL_HH

#define KURAMOTOMODEL_HH
#include "sol_rk4d.hh"
#define _USE_MATH_DEFINES
#include <math.h>
#include <random> 
#include <omp.h>
#include <cstring>
#include <cstdio>

class KuramotoModel_Solver : public rk4d{
// class principle:
// q[0:N-1] contains x states
// q[N:2N-1] contains y states
// q[2N: 3N-1] contains phase 
	public:
		int N;
		double* vi;
		double* wi; 
		double J;
		double K;
		double A; 
		double B;
		double* timeShows;
		int len_timeShows;
		bool printSelect;
		
		
	/* Construct the KuromotoModel */
		KuramotoModel_Solver(
			int N_, double* timeShow, int len_timeShow, bool printSel, double J_, double K_) : 
			N(N_), vi(new double[N_*2]), wi(new double[N_]), J(J_), K(K_),
			timeShows(timeShow), len_timeShows(len_timeShow), printSelect(printSel),
			rk4d(3*N_){};


	/* Define Initialization function for ODE */
		virtual void init(){
			//Initialize constants A = B = 1 ; vi = wi = 0
			A=B=1;
			for(int i=0; i<N; i++){
				vi[i] = i/(1250.);// originally set to 0
				vi[i+N] = i/(1000.); // originally set to 0
				wi[i] = i/500.; // originally set to 0
			}
			
			// Initialize (q) xi,yi as random in unit desk and
			// the phases as random 0 to 2pi
			static std::default_random_engine e;
			static std::uniform_real_distribution<> dis(0, 1);
			for(int i=0; i<N; i++){
				double sig = (double)(dis(e)>0.5);
				double x = (double)dis(e);
				x = sig*x - (1.-sig)*x;
				sig = (double)(dis(e)>0.5);
				double y = (double)dis(e)*sqrt(1-x*x);
				y = sig*y - (1.-sig)*y;

				q[i] = x;
				q[i+N] = y;
				q[i+2*N] = double(dis(e))*2*M_PI; 
			}
		};
	
	
	/* Function to print the grid */
		virtual void print(double t_, double *in, double error){
			// In = state vector q
			// printSelect true causes only len_timeShow number of prints to execute
			// at the timestamps where t_ = timeShows
			// Assumed, driver code calls printDense true and 
			// that timeshows aligns with dense step
			if(printSelect){
				double dostatus = 0;
				for(int i=0;i<len_timeShows;i++){
					if(t_==timeShows[i]){
						dostatus = 1;
						break;
					}
				}
				if(dostatus==0) return;
			}

			for(int i=0; i<N; i++){
				printf("%g %g %g ", in[i], in[i+N], in[i+2*N]);
				puts("");
			}
			puts("\n");
			
		};   
	
	
	/* define derivative function */
		virtual void ff(double t_,double *in,double *out){
			double* dqi = new double[3*N];
			#pragma omp parallel for
			for(int i = 0; i<N; i++){
				double* dxhold = new double[N];
				double* dyhold = new double[N];
				double* dthetahold = new double[N];
				#pragma omp parallel for
				for(int j=0; j<N; j++){
					if(i==j){
						dxhold[j] = 0.;
						dyhold[j] = 0.;
						dthetahold[j] = 0.;
					}else{
						double dx = in[j]-in[i];
						double dy = in[j+N]-in[i+N];
						double dtheta = in[j+2*N]-in[i+2*N];
						double scale = (A+J*cos(dtheta));
						double mag = (dx*dx+dy*dy);
						
						dxhold[j] = dx/sqrt(mag)*scale - B*dx/(mag);
						dyhold[j] = dy/sqrt(mag)*scale - B*dy/(mag);
						dthetahold[j] = sin(dtheta)/mag;
					}
				}
				// can parallelize this addition but dont need to
				double dxsum = 0.;
				double dysum = 0.;
				double dthetasum = 0.;
				for(int k=0; k<N; k++){
					dxsum+=dxhold[k];
					dysum+=dyhold[k];
					dthetasum+=dthetahold[k];
				}
				dqi[i] = vi[i] + dxsum/N;
				dqi[i+N] = vi[i+N] + dysum/N;
				dqi[i+2*N] = wi[i] + K*dthetasum/N; 
				delete [] dxhold;
				delete [] dyhold;
				delete [] dthetahold;
			}
			memcpy(out,dqi,3*N*sizeof(double));
			delete [] dqi;
		};
	
};
#endif



















