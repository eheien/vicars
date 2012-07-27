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
  const double log_min;
  const double poly_pow;
  const double poly_c;
  const int npoints_poly;
  const int npoints_log;
  const double log_takeover;
  double log_dip_polynomial(double v) const { 
    return (-1 * log_min * poly_c * pow(v,poly_pow) + log_min);
  };
  double intersection(double a, double b);

public:

  double operator ()(double v) const {
    if (abs(v) < log_takeover/2) {
    return gsl_spline_eval(spline,v,acc); 
    } else { return log(abs(v)); }
  };

  double deriv(double v) const {
    if (abs(v) < log_takeover/2) {
    return gsl_spline_eval_deriv(spline,v,acc); 
    } else { return abs(1/v); };
  };
    
  LogSpline(double _log_min, double _poly_pow, double _poly_c,
      int _npoints_poly, int _npoints_log, double _log_takeover);
  ~LogSpline(void);
};

