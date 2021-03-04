// Define the Brusselator problem
#ifndef BRUSSELATOR_HH
#define BRUSSELATOR_HH
#include "sol_rk4d.hh"
#include "sol_rk4classic.hh"
#include <cstdio>

/** This class has functions to specify the test Brusselator problem. */
class brusselator {
    public:
        /** Evaluates the function f(x,y) on the RHS of the ODE.
         * \param[in] t the dependent x variable in Hairer et al.'s notation.
         * \param[in] in the array containing y variable.
         * \param[in] out the function f(x,y). */
        inline void brus_ff(double t_,double *in,double *out) {
            double &y1=*in,&y2=in[1];
            *out=1+y1*(y1*y2-4);
            out[1]=y1*(3-y1*y2);
        }
        /** Sets up the initial conditions for the ODE.
         * \param[in] q the array to write to. */
        inline void brus_init(double *q) {
            *q=1.5;
            q[1]=3.;
        }      
};

/** Class to solve the Brusselator problem with the fourth-order Runge-Kutta method. */
class brus_rk4 : public rk4d, public brusselator {
    public:
        brus_rk4() : rk4d(2) {}
        virtual void ff(double t_,double *in,double *out) {
            brus_ff(t_,in,out);
        }
        virtual void init() {brus_init(q);}
        virtual void print(double t_,double *in, double error){
				printf("%g", t_);
				for(int i=0;i<dof;i++)printf(" %g",in[i]);
				printf(" %g", error);
				puts("");
		  }
};

/** compare to the classic RK4 method */
class brus_rk4_reference : public classic_rk4, public brusselator {
    public:
        brus_rk4_reference() : classic_rk4(2) {}
        virtual void ff(double t_,double *in,double *out) {
            brus_ff(t_,in,out);
        }
        virtual void init() {brus_init(q);}
        virtual void print(double t_,double *in){
				printf("%g", t_);
				for(int i=0;i<dof;i++)printf(" %g",in[i]);
				puts("");
		  }
};


#endif
