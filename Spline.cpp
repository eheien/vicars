#include "Spline.h"

LogSpline::LogSpline(void) : log_min(LOGMIN) {

  unsigned int i;
  std::vector<double> v, logv;
  double vi, pv, pweight, step;

  double isct = intersection(log_dip_polynomial, log, 1e-12,1);

  v.push_back(0);
  logv.push_back(log_min);

  pweight = 1;
  step = isct/10;
  for (vi=step; vi < isct; vi+=step) {
    v.push_back(vi);
    pv = log_dip_polynomial(vi);
    logv.push_back(pweight * pv + (1-pweight) * log(vi));
    pweight *= 0.75;
  }

  pweight = 1.01;
  while (vi < 1e12) {
    v.push_back(vi);
    logv.push_back(log(vi));
    vi *= pweight;
    pweight *= 1.02;
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

  acc = gsl_interp_accel_alloc ();
  spline = gsl_spline_alloc (gsl_interp_akima, num_points);
  gsl_spline_init (spline, v_arr, logv_arr, num_points);

#ifdef OUTPUT
  std::ofstream log_approx_graph;
  log_approx_graph.open ("log_approx.txt");
  for (vi=10; vi<1e11; vi *= 10) {
    log_approx_graph << vi << " " << gsl_spline_eval(spline, vi, acc)  << std::endl;
  }
  log_approx_graph.close();
#endif

#ifdef DEBUG
  std::cout << "POW/LOG INTERSECTION : " << isct << std::endl;
  std::cout << "# SAMPLE POINTS : " << v.size() << std::endl;
#endif
}

LogSpline::~LogSpline(void) {
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
}

double intersection(double (*f)(double), double (*g)(double), double a, double b) {
  unsigned int i = 1;
  double fgp,fga,p;
  fga = f(a) - g(a);

  while (i<MAXITER) {
    p = a + (b-a)/2;
    fgp = f(p) - g(p);
    if (fgp == 0 or (b-a)/2 < TOL) break;
    i++;
    if (fga*fgp > 0) {
      a = p;
      fga = fgp;
    } else { b = p; };
  }
  return p;
};

double log_dip_polynomial(double v) { return (-1 * LOGMIN * 1e13 * pow(v,2.5) + LOGMIN); };
