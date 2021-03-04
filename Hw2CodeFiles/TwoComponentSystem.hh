// Define the Brusselator problem
#ifndef TWOCOMPONENTSYSTEM_HH
#define TWOCOMPONENTSYSTEM_HH
#include "sol_rk4d.hh"
#include "sol_IRK5.hh"
#include <cstdio>
#include <cmath>

// This class is created to specifically test the dense output functionality
// according to problem 1B in the homework
class TwoComponentSystem{
    public:
        inline void this_ff(double t_,double *in,double *out) {
            double &y1=*in,&y2=in[1];
            *out=-1.*t_*y2;
            out[1]=1.*t_*y1;
        }
        
        /** Sets up the initial conditions for the ODE.
         * \param[in] q the array to write to. */
        inline void this_init(double *q) {
            *q=1.;
            q[1]=0.;
        }
};


/** Class to solve the two-component system with the fourth-order Runge-Kutta */
class twoComponent_rk4 : public rk4d, public TwoComponentSystem {
    public:
        twoComponent_rk4() : rk4d(2) {}
        virtual void ff(double t_,double *in,double *out){this_ff(t_,in,out);}
        virtual void init(){this_init(q);}
        virtual void print(double t_,double *in, double error){
				printf("%g", t_);
				for(int i=0;i<dof;i++)printf(" %g",in[i]);
				// this special print adds the errors to true sol
				printf(" %g", in[0]-cos(t_*t_/2.) );
				printf(" %g", in[1]-sin(t_*t_/2.) );
				puts("");
			}
};

/** Class to solve the two-component system with the fifth order IRK5 */
class twoComponent_irk5 : public irk5, public TwoComponentSystem{
	public: 
		twoComponent_irk5() : irk5(2){}
		virtual void ff(double t_,double *in,double *out) {
            this_ff(t_,in,out);
        }
		virtual void init() {this_init(q);}
		virtual void print(double t_,double *in){
				printf("%g", t_);
				for(int i=0;i<dof;i++)printf(" %g",in[i]);
				// this special print adds the errors to true sol
				printf(" %g", in[0]-cos(t_*t_/2.) );
				printf(" %g", in[1]-sin(t_*t_/2.) );
				puts("");
		}
};


#endif
