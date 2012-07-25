#include "Spline.h"
#include <iostream>

LogSpline::LogSpline(void) {

  unsigned int i;
  const double log_min = -20;
  std::vector<double> v, logv;
  double vi;

  v.push_back(0);
  logv.push_back(log_min);

  for (vi=exp(log_min) + 1e-10; vi < 1e5; vi *= 1.8) {
    v.push_back(vi);
    logv.push_back(log(vi));
  }

  std::vector<double> vrev = v;
  reverse(vrev.begin(),vrev.end());
  vrev.pop_back();
  for (i=0; i<vrev.size(); ++i) {
    vrev[i] = -vrev[i];
  }

  std::vector<double> logvrev = logv;
  reverse(logvrev.begin(),logvrev.end());
  logvrev.pop_back();

  v.insert(v.begin(), vrev.begin(), vrev.end());
  logv.insert(logv.begin(), logvrev.begin(), logvrev.end());

  size_t num_points = v.size();
  double* v_arr = &v[0];
  double* logv_arr = &logv[0];
	
	//for (int i=0;i<num_points;++i) std::cerr << v_arr[i] << " " << logv_arr[i] << std::endl;

  acc = gsl_interp_accel_alloc ();
  spline = gsl_spline_alloc (gsl_interp_cspline, num_points);
  gsl_spline_init (spline, v_arr, logv_arr, num_points);

	/*std::cerr << std::endl;
	std::cerr << std::endl;
	for (vi=exp(log_min) + 1e-10; vi < 1e5; vi *= 1.8) {
		std::cerr << vi << " " << gsl_spline_eval(spline,vi,acc) << std::endl;
	}*/
}

LogSpline::~LogSpline(void) {
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
}
