#ifndef STARORBIT_HH
#define STARORBIT_HH

#include "sol_IRK5.hh"
#include <math.h>
#include <cstdio>

// Simulate the dynamics of a star's movement in the galaxy
class orbit_irk5 : public irk5{
	public:
	/* Initialize the star galaxy orbit problem */
		orbit_irk5() : irk5(6){}
		
	/* set the derivative function */
		virtual void ff(double t_,double *in,double *out){ 
			// constants for our particular problem
			double omega=0.25, A=1., C=1., a=1.25, b=1., c=0.75;

			// computes the derivatives for the position and momentum
			// after dt time step
			// in is state vector q holding 6 dof: p1, p2, p3, q1, q2, q3
			// out is the change of the quantity computed by the 
			// derivative of the hamiltonian
			double &p1 = in[0], &p2 = in[1], &p3 = in[2];
			double &q1 = in[3], &q2 = in[4], &q3 = in[5];
			double phi = C + q1*q1/(a*a) + q2*q2/(b*b)+ q3*q3/(c*c);
			
			out[0] = -1.*(-omega*p2 + A/phi*(2.*q1)/(a*a));
			out[1] = -1.*(omega*p1 + A/phi*(2.*q2)/(b*b));
			out[2] = -1.*(A/phi*(2.*q3)/(c*c));
			out[3] = p1+omega*q2;
			out[4] = p2-omega*q1;
			out[5] = p3;
		}
		
	/* set the initial conditions */
		virtual void init(){
			// constants for our particular problem
			double omega=0.25, A=1., C=1., a=1.25, b=1., c=0.75;
			
			// initialize the starting state vector contained in q
			// q = [p1, p2, p3, q1, q2, q3]
			// set p1(0) = 0, p3(0) = 0.2
			// set q1(0) = 2.5, q2(0)=0, q3(0)=0
			// p2 is obtained by solving root and taking larger value
			q[0] = 0.;
			q[2] = 0.2;
			q[3] = 2.5;
			q[4] = 0.;
		 	q[5] = 0.;
			double Hzero = 2.;
			
			// solve for p2 with H(0)=2
			double p1 = q[0], p3 = q[2];
			double q1 = q[3], q2 = q[4], q3 = q[5];		
			double z0 = A*log(C + q1*q1/(a*a)+ q2*q2/(b*b) + q3*q3/(c*c) ) + omega*p1*q2 + p1*p1/2. + p3*p3/2. - Hzero;
			double z1 = -omega*q1;
			double z2 = 1./2.;
			double r1 = (-z1 + sqrt(z1*z1 - 4.*z2*z0))/(2.*z2);
			double r2 = (-z1 - sqrt(z1*z1 - 4.*z2*z0))/(2.*z2);
			//printf("Root for Initialize P2: %g %g\n", r1, r2);
			if(r1>r2) q[1]=r1;
			if(r2>r1) q[1]=r2;
		}
		
	/* compute the hamiltonian */
		inline double hamiltonian(){
			double omega=0.25, A=1., C=1., a=1.25, b=1., c=0.75;
			double p1 = q[0], p2 = q[1], p3 = q[2];
			double q1 = q[3], q2 = q[4], q3 = q[5];
			return (p1*p1 + p2*p2 + p3*p3)/2. + omega*(p1*q2-q1*p2) + A*log(C + q1*q1/(a*a)+ q2*q2/(b*b) + q3*q3/(c*c) );

		}
		
   /* set the print function */
		virtual void print(double t_,double *in){ 
			printf("%g %g ", t_, hamiltonian());
			for(int i=0; i<dof; i++) printf("%g ", in[i]);
			puts("");
		}
};

#endif























