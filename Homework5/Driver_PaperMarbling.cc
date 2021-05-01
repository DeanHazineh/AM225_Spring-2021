#include <sys/types.h>
#include <sys/stat.h>
#include <cmath>
#include "fluid_2d.hh"

void runSolver(fluid_2d &fd, double t, int frames){
	fd.set_boundaries();    
	fd.initialize(0,0.6);	 
	fd.solve(t, frames);	
}

int main() {
// View vector fields at 
//https://www.geogebra.org/m/cXgNb58T


	/* Create simulation class */
    const char fn[]="paperMarble.out";
    mkdir(fn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
    unsigned int fflags=0; //1: horizontal velocity, 2:vertical velocity, 4: pressure
    const double ax=-1., bx=1.;
    const double ay=-1., by=1.;
    fluid_2d f2d(256,256,false,false,ax,bx,ay,by,0.002,1.,fflags,fn);
	 f2d.resetVelocityPressureField();

	  
	// First Layer - Specify an initial color and field distribution the solve
	 f2d.init_fields_CustomFunc();
	 f2d.init_Color(); // creates black and purple grid
	 const char fn1[]="paperMarble.out/colorField1.png";
	 f2d.outputColorField(fn1);
	 runSolver(f2d, 3.5, 5);


	// Second Layer -- use same advection field as before
	 int numCircles=35;
	 double maxRad = 0.1, useRad=0., usex=0., usey=0.;
	 for(int i=0;i<numCircles;i++){
	 	useRad = maxRad/RAND_MAX*double(rand());
	 	usex = ax+(bx-ax)/RAND_MAX*double(rand());
	 	usey = ay+(by-ay)/RAND_MAX*double(rand());
	 	 f2d.AddCircle_colorfield(usex, usey, 1.2*useRad, 34, 131, 164);
	 	 f2d.AddCircle_colorfield(usex, usey, 0.8*useRad, 19, 55, 105);
	 } 
	 const char fn2[]="paperMarble.out/colorField2.png";
	 f2d.outputColorField(fn2);
 	 runSolver(f2d, 3.5, 5);
	
	
	// Third Layer -- add color and new field
	 f2d.resetVelocityPressureField();
	 f2d.renew_VelocityField1(0, 0, 0.3);
	 f2d.AddCircle_colorfield(0., 0., 0.35, 248, 188, 4);
	 f2d.AddCircle_colorfield(0., 0.1, 0.1, 164, 67, 34);
    f2d.AddCircle_colorfield(0., -0.1, 0.1, 164, 67, 34);
    f2d.AddCircle_colorfield(0., 0.0, 0.2, 164, 67, 34);
    const char fn3[]="paperMarble.out/colorField3.png";
    f2d.outputColorField(fn3);
  	 runSolver(f2d, 3.5, 5);
	
	// Fourth Layer -- Similar to above but a different advection current
	 f2d.resetVelocityPressureField();
    f2d.AddCircle_colorfield(0., 0.0, 0.35, 94, 2, 124);
    f2d.AddCircle_colorfield(0., 0.0, 0.2, 188, 4, 248);
    f2d.init_fields();
    f2d.renew_VelocityField2(0, 0, 0.3);
	 const char fn4[]="paperMarble.out/colorField4.png";
    f2d.outputColorField(fn4);
  	 runSolver(f2d, 3.5, 5);
	
	// Fifth Layer -- 
	 numCircles=50;
	 maxRad = 0.05, useRad=0., usex=0., usey=0.;
	 for(int i=0;i<numCircles;i++){
	 	useRad = maxRad/RAND_MAX*double(rand());
	 	usex = ax+(bx-ax)/RAND_MAX*double(rand());
	 	usey = ay+(by-ay)/RAND_MAX*double(rand());
	 	 f2d.AddCircle_colorfield(usex, usey, 1.2*useRad, 201, 252, 244);
	 	 f2d.AddCircle_colorfield(usex, usey, 0.8*useRad, 241, 254, 252);
	 } 
	 const char fn5[]="paperMarble.out/colorField5.png";
	 f2d.outputColorField(fn5);
 	 runSolver(f2d, 0.5, 3);
	
	

}
