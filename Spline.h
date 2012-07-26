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
#define OUTPUT

class LogSpline {
private:		
  gsl_interp_accel *acc;
  gsl_spline *spline;
  const double log_min;
  const double pow;

public:
  double operator ()(double v) { return gsl_spline_eval(spline,v,acc); };
  double deriv(double v) { return gsl_spline_eval_deriv(spline,v,acc); };
  LogSpline(void);
  ~LogSpline(void);
};

double intersection(double (*f)(double), double (*g)(double), double a, double b);
double log_dip_polynomial(double v);
