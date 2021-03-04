#ifndef DISCONTINUOUSTWOCOMPONENT_HH
#define DISCONTINUOUSTWOCOMPONENT_HH

#include "sol_rk4classic.hh"
#include "sol_rk4d.hh"
#include <cstdio>
#include <cmath>

/* Create a class that encompasses the ODE in Problem 7 */
class discontinuous_twocomponent{
    public:
		  /* Formula to update Derivatives for the system */
        inline void this_ff(double t_,double *in, double *out) {
        		// in contains the state x, y and out stores the derivative
        		double &x=in[0], &y=in[1];
        		if( fabs(x) >= fabs(y) ){
        			*out = 0.;
        			out[1] = x;
        		}else if(fabs(y)>fabs(x)){
        			*out = -y;
        			out[1] = 0.;
        		}
        }
        
        /** Sets up the initial conditions for the ODE.*/
        inline void this_init(double *q) {
            *q=1.;
            q[1]=0.;
        }
        
        /* define analytic solution x */
        inline double xanalytic(double t){
				//period is exactly 8 so define effective time in periodic signal
				double teff = fmod(t,8.);
				if(teff>=0. && teff<=1.) return 1.;	
				else if(teff>1. && teff<3.) return 1.-(teff-1.);
				else if(teff>=3. && teff<=5.) return -1.;
				else if(teff>5. && teff<7.) return -1. + (teff-5.);
				else if(teff>=7. && teff<=9.) return 1.;
				else{
					printf("ERROR in Analytic Computation X");
					exit(EXIT_FAILURE);
				}
			}
			
		/* define analytic soultion y */
		inline double yanalytic(double t){
			//period is exactly 8 so define effective time in periodic signal
			double teff = fmod(t,8.);
			// define the analytic step equations
			if(teff>=0. && teff <=1.) return 0.+teff;	
			else if(teff>1. && teff<3.) return 1.;
			else if(teff>=3. && teff<=5.) return 1.-(teff-3.);
			else if(teff>5. && teff<7.) return -1.;
			else if(teff>=7. && teff<=9.) return -1. + (teff-7.);
			else{
				printf("ERROR in Analytic Computation Y");
				exit(EXIT_FAILURE);
			}
		}
};


/* Derived solver class which steps the ODE using classic rk4 Algorithm*/
class discontinuous_twocomponent_rk4classic :  public classic_rk4, public discontinuous_twocomponent{
	public:
     discontinuous_twocomponent_rk4classic() : classic_rk4(2) {}

     virtual void init(){ this_init(q); }
     virtual void ff(double t_,double *in,double *out){ this_ff(t_, in, out); }   
     virtual void print(double t_,double *in){
		  	// print time of state
		  	printf("%g ", t_);
		  	// print state value from ode solver
		  	for(int i = 0; i<dof; i++) printf("%g ", in[i]);
		  	// print the analytic state value
		  	printf("%g %g ", xanalytic(t_), yanalytic(t_));
		  	puts("");
     }
};

/* derived solver class which steps the ODE using FSAL RK4(3) dynamic step size */
class discontinuous_twocomponent_rk4d : public rk4d, public discontinuous_twocomponent{
	public:
     discontinuous_twocomponent_rk4d() : rk4d(2) {}

     virtual void init(){ this_init(q); }
     virtual void ff(double t_,double *in,double *out){ this_ff(t_, in, out); }   
     virtual void print(double t_,double *in, double error){
		  	// print time of state
		  	printf("%g ", t_);
		  	// print state value from ode solver
		  	for(int i = 0; i<dof; i++) printf("%g ", in[i]);
		  	// print the analytic state value
		  	printf("%g %g ", xanalytic(t_), yanalytic(t_));
		  	puts("");
     }
};

#endif
