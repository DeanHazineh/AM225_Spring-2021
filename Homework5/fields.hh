/*Edited by Dean for AP 225 Assignment 5
Introduced C vector R,G,B as part of the structure */ 

#ifndef FIELDS_HH
#define FIELDS_HH

#include <cmath>

/** Data structure for storing the fields at grid points. */
struct field {
	 /** Introduce Color values as part of the field distribution  */
	 double R;
	 double G;
	 double B;

    /** The horizontal velocity. */
    double u;
    /** The vertical velocity. */
    double v;
    /** The pressure. */
    double p;
    /** The intermediate horizontal velocity. */
    double us;
    /** The intermediate vertical velocity. */
    double vs;
    
    inline void prd_bc(field &f) {
        u=f.u;v=f.v;
    }
    inline void no_slip(field &f) {
        u=-f.u;v=-f.v;
    }
    inline void no_slipColor(field &f){
    	R = f.R; G= f.G; B = f.B;
    }
    
    /** Computes the maximum allowable timestep based on the CFL restriction
     * from the velocity stored in this class.
     * \param[in] (xsp,ysp) the horizontal and vertical grid spacings.
     * \return The reciprocal of the maximum allowable timestep. */
    inline double cfl(double xsp,double ysp) {
        double uc=fabs(u)*xsp,vc=fabs(v)*ysp;
        return uc>vc?uc:vc;
    }
};

#endif
