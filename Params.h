#include <math.h>

#ifndef _VICARS_PARAMS_H_
#define _VICARS_PARAMS_H_

class SimParams {
public:
    double                  V_star; //
	double                  G;		// shear modulus parameter, pascals
	double                  D_c;    // characteristic slip, meters
	double                  W;
	double                  beta;   // shear wave speed, meters/second
    double                  sigma;  // effective normal stress, pascals
    double                  A;
    double                  B;
    double                  v_ss;   // steady state pulling speed, meters/second
    double                  mu_0;   // nominal coefficient of friction
    double                  side;   // side length, meters
    double                  rho;   // density of material (kg/m^3)
    double                  gravity;   // acceleration due to gravity (m/s^2)
    double                  L;      // rate state length normalization
    
    double area(void) const { return side*side; };
    double mass(void) const { return rho*pow(side, 3); };
    double k(void) const { return G*sqrt(area())/2; };
    double tau(void) const { return 2*M_PI*sqrt(mass()/k()); };
};

#endif
