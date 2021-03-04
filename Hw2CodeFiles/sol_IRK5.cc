
#include "sol_IRK5.hh"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

// These function calls are inspired by the sol_hammer_h example
// provided by Rycroft in the 1a_ode_solvers package since it is a similar implicit
// method

irk5::irk5(int dof_) : dof(dof_), fcount(0), t(0.),
	 q(new double[dof]),	dq(new double[dof]),
    k1(new double[dof]), k2(new double[dof]), k3(new double[dof]),
    k1b(new double[dof]), k2b(new double[dof]), k3b(new double[dof])
    {
		// Initializes the fourth-order Runge-Kutta solver
		// with fixed parameterization for facx
    }

irk5::~irk5() {
    delete [] k3b;
    delete [] k2b;
    delete [] k1b;
    delete [] k3;
    delete [] k2;
    delete [] k1;
    delete [] dq;
    delete [] q;
}

bool irk5::solve_Poincaire(double duration, int steps){
    // Set up initial condition and compute timestep    
    init();
    double dt=duration/steps;
    // hold previous step in a special array
	 double* prevStep = new double[dof];
 	 memcpy(prevStep,q,dof*sizeof(double));
	
    // Perform integration steps
    for(int i=0;i<steps;i++) {
        if(!step(dt)) {
            fputs("Too many iterations in IRK method\n",stderr);
            return false;
        }
        // check poincaire condition
		  if(prevStep[4]<0 && q[4]>0){//zero crossing in y (because want xz plane)
		  		if(q[3]>0){// in +x plane
		  			// do simple linear interpolation for crossing
		  			double xlin = prevStep[3] + (q[3]-prevStep[3])/(q[4]-prevStep[4])*(0.-prevStep[4]);
		  			double zlin = prevStep[5] + (q[5]-prevStep[5])/(q[4]-prevStep[4])*(0.-prevStep[4]);
		  			printf("%g %g %g",t, xlin, zlin); puts("");
	  			}	 	
		  }   
		  memcpy(prevStep,q,dof*sizeof(double));   
    }
    return true;
}

bool irk5::solve_fixed(double duration, int steps, bool output){
    // Set up initial condition and compute timestep
    init();
    double dt=duration/steps;

    // Perform integration steps
    if(output) print(t, q);
    for(int i=0;i<steps;i++) {
        if(!step(dt)) {
            if(output) fputs("Too many iterations in IRK method\n",stderr);
            return false;
        }
        if(output) print(t,q);
    }
    return true;
}

bool irk5::step(double dt) {
    int iter=0;
    double delsq, d,*c;
    // precompute the coefficients required for solving
    const double r1m = (16.-sqrt(6))/72.;
    const double r1p = (16.+sqrt(6))/72.;
    const double r2m = (328.-167.*sqrt(6))/1800.;
    const double r2p = (328.+167.*sqrt(6))/1800.;
    const double r3m = (-2.-3.*sqrt(6))/450.;
    const double r3p = (-2.+3.*sqrt(6))/450.;
    const double r4m = (85.-10.*sqrt(6))/180.;
    const double r4p = (85.+10.*sqrt(6))/180.;
    const double tm = (4.-sqrt(6))/10.;
    const double tp = (4.+sqrt(6))/10.;
    const double c1m = (16.-sqrt(6))/36.;
    const double c1p = (16.+sqrt(6))/36.;

    // Clear steps
    for(int i=0;i<dof;i++) k1[i]=k2[i]=k3[i]=0;
    do {
        // Check for too many iterations
        if(++iter>1000) {
            return false;
            //fputs("Too many iterations in IRK\n",stderr);
            //exit(1);
        }

        // Perform update
        for(int i=0;i<dof;i++) dq[i]=q[i] + dt*(r1m*k1[i] + r2m*k2[i] + r3p*k3[i]);
        ff(t+dt*tm,dq,k1b);
        
        for(int i=0;i<dof;i++)  dq[i]=q[i] + dt*(r2p*k1[i] + r1p*k2[i] + r3m*k3[i]);
        ff(t+dt*tp,dq,k2b);
      
        for(int i=0;i<dof;i++)  dq[i]=q[i] + dt*(r4m*k1[i] + r4p*k2[i] + (1./18.)*k3[i]);
        ff(t+dt*1.,dq,k3b);
        
        fcount+=3;

        // Find size of step from previous iteration
        delsq=0;
        for(int i=0;i<dof;i++) {
            d=k1[i]-k1b[i];delsq+=d*d;
            d=k2[i]-k2b[i];delsq+=d*d;
            d=k3[i]-k3b[i];delsq+=d*d;
        }

        // Switch ki<->kib  array pointers. This will make ki
        // be used as the new values on the next iteration.
        c=k1b;k1b=k1;k1=c;
        c=k2b;k2b=k2;k2=c;
        c=k3b;k3b=k3;k3=c;
    } while(delsq>1e-25);

    // Complete solution
    for(int i=0;i<dof;i++) q[i]+=dt*(c1m*k1[i] + c1p*k2[i] + (1./9.)*k3[i]);
    t+=dt;
    return true;
}

