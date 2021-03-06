#include "Spline.h"
#include <map>

BasicSpline::BasicSpline(double trans_value) : transition_value(trans_value) {
	std::map<double, double> v;
	std::map<double, double>::const_iterator	it;
	unsigned int	i, num_steps = 20;
	double			step, val, lower_bound = 1e-10;
	
	step = pow((trans_value/lower_bound),1.0/num_steps);
	for (val=trans_value,i=0;i<num_steps;i++,val/=step) {
		v[val] = v[-val] = log(val);
	}
	v[0] = log(lower_bound);
	
	size_t num_points = v.size();
	double *v_arr = new double[num_points], *logv_arr = new double[num_points];
	for (i=0,it=v.begin();it!=v.end();++it,++i) {
		v_arr[i] = it->first;
		logv_arr[i] = it->second;
	}
	
	acc = gsl_interp_accel_alloc ();
	spline = gsl_spline_alloc (gsl_interp_akima, num_points);
	gsl_spline_init (spline, v_arr, logv_arr, num_points);
	
#ifdef OUTPUT_SPLINE
	std::ofstream log_sample_points;
	unsigned int sample_steps = 250;
	step = pow((trans_value/lower_bound),1.0/sample_steps);
	log_sample_points.open ("log_sample_points.txt");
	for (val=trans_value,i=0;i<sample_steps;i++,val/=step) {
		log_sample_points << val << " " << gsl_spline_eval(spline, val, acc) << " " << log(fabs(val)) << std::endl;
	}
	log_sample_points.close();
#endif
	delete[] v_arr;
	delete[] logv_arr;
}

BasicSpline::~BasicSpline(void) {
	gsl_spline_free (spline);
	gsl_interp_accel_free (acc);
}

LogSpline::LogSpline(double _log_min, double _poly_pow, double _poly_c, double _log_takeover, int _npoints_poly) : log_min(_log_min), poly_pow(_poly_pow), poly_c(_poly_c), log_takeover(_log_takeover), npoints_poly(_npoints_poly) {
	
	unsigned int i;
	std::vector<double> v, logv;
	double vi, pv, pweight, step;
	
	double isct = intersection(1e-12,1);
	
	v.push_back(0);
	logv.push_back(log_min);
	
	pweight = 1;
	step = isct/npoints_poly;
	for (vi=step; vi < isct; vi+=step) {
		v.push_back(vi);
		pv = log_dip_polynomial(vi);
		logv.push_back(pweight * pv + (1-pweight) * log(vi));
		pweight *= 0.5;
	}
	
	pweight = 1.01;
	while (vi < log_takeover) {
		v.push_back(vi);
		logv.push_back(log(vi));
		pweight *= 1.03;
		vi *= pweight;
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
	spline = gsl_spline_alloc (gsl_interp_cspline, num_points);
	gsl_spline_init (spline, v_arr, logv_arr, num_points);
	
#ifdef OUTPUT_SPLINE
	std::ofstream log_approx_graph, log_sample_points;
	log_approx_graph.open ("log_approx.txt");
	log_sample_points.open ("log_sample_points.txt");
	for (vi=-1e-3; vi<1e-3; vi += (1e-3)/1000) {
		log_approx_graph << vi << " " << gsl_spline_eval(spline, vi, acc)  << std::endl;
	}
	for (i=0; i<v.size(); i++) {
		log_sample_points << v[i] << " " << logv[i] << std::endl;
	}
	log_approx_graph.close();
	log_sample_points.close();
#endif
	
#ifdef DEBUG_SPLINE
	std::cout << "POW/LOG INTERSECTION : " << isct << std::endl;
	std::cout << "# SAMPLE POINTS : " << v.size() << std::endl;
#endif
}

LogSpline::~LogSpline(void) {
	gsl_spline_free (spline);
	gsl_interp_accel_free (acc);
}

double LogSpline::intersection(double a, double b) {
	unsigned int i = 1;
	double fgp,fga,p;
	fga = log_dip_polynomial(a) - log(a);
	
	while (i<MAXITER) {
		p = a + (b-a)/2;
		fgp = log_dip_polynomial(p) - log(p);
		if (fgp == 0 or (b-a)/2 < TOL) break;
		i++;
		if (fga*fgp > 0) {
			a = p;
			fga = fgp;
		} else { b = p; };
	}
	return p;
};
