#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <algorithm>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>

#define TOL 1e-15
#define MAXITER 1000

//#define DEBUG
//#define OUTPUT

class LogSpline {
private:		
  gsl_interp_accel *acc;
  gsl_spline *spline;
  const double log_min, poly_pow, poly_c, log_takeover;
  const int npoints_poly;

  double log_dip_polynomial(double v) const { 
    return (-1 * log_min * poly_c * pow(v,poly_pow) + log_min);
  };

  double intersection(double a, double b);

public:

  double operator ()(double v) const {
    if (fabs(v) < 0.9*log_takeover) {
    return gsl_spline_eval(spline,v,acc); 
    } else { return log(fabs(v)); }
  };

  double deriv(double v) const {
    if (fabs(v) < 0.9*log_takeover) {
    return gsl_spline_eval_deriv(spline,v,acc); 
    } else { return fabs(1/v); };
  };
    
  LogSpline(double _log_min, double _poly_pow, double _poly_c,
            double _log_takeover, int _npoints_log);
  ~LogSpline(void);
};

