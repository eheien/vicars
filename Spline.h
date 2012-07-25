#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <algorithm>
#include <math.h>
#include <vector>

class LogSpline {
private:		
  gsl_interp_accel *acc;
  gsl_spline *spline;
  double log_min;

public:
  double operator ()(double v) { return gsl_spline_eval(spline,v,acc); };
  LogSpline(void);
  ~LogSpline(void);
};
