#include "sol_rk4classic.hh"
#include <cstdio>
#include <cstring>
#include <math.h>  

classic_rk4::classic_rk4(int dof_) : dof(dof_), fcount(0), t(0.), q(new double[dof]),
    dq(new double[dof]),
    k1(new double[dof]), k2(new double[dof]), k3(new double[dof]),
    k4(new double[dof]) {}

/** The class destructor frees the dynamically allocated memory. */
classic_rk4::~classic_rk4() {
    delete [] k4;
    delete [] k3;
    delete [] k2;
    delete [] k1;
    delete [] dq;
    delete [] q;
}

/** Performs an integration step with the fourth-order Runge-Kutta solver.
 * \param[in] dt the integration step. */
void classic_rk4::step(double dt) {

    // Second RK step
    for(int i=0;i<dof;i++) dq[i]=q[i]+0.5*dt*k1[i];
    ff(t+0.5*dt,dq,k2);

    // Third RK step
    for(int i=0;i<dof;i++) dq[i]=q[i]+0.5*dt*k2[i];
    ff(t+0.5*dt,dq,k3);

    // Fourth RK step
    t+=dt;fcount+=4;
    for(int i=0;i<dof;i++) dq[i]=q[i]+dt*k3[i];
    ff(t,dq,k4);

    // Complete solution
    for(int i=0;i<dof;i++) dq[i]=q[i]+dt*(1/6.)*(k1[i]+2*(k2[i]+k3[i])+k4[i]);

    // Reuse k2 to store the derivative at the new solution
    ff(t,dq,k2);
}


/** Solves the ODE problem with a fixed integration step.
 * \param[in] duration the time to solve for.
 * \param[in] steps the number of integration steps
 * \param[in] output whether to print each integration step. */
void classic_rk4::solve_fixed(double duration,int steps,bool output) {

    init();
    double dt=duration/steps;

    // Perform integration steps
    if(output) print(t,q);
    ff(t,q,k1);
    for(int i=0;i<steps;i++) {
        step(dt);

        // Move the required data for the next step into position
        memcpy(q,dq,dof*sizeof(double));
        memcpy(k1,k2,dof*sizeof(double));
        if(output) print(t,q);
    }
}
