
#include "RateState.h"

#define NBLOCKS			1
//#define USE_SIMPLIFIED_EQNS

int main(int argc, char **argv)
{
	ViCaRS			sim(NBLOCKS);
	unsigned int	i;
	double			param_a, param_b, param_k, param_r;
	FILE			*fp;
	int				res;
	
	realtype		x_0, v_0, h_0, x_err, v_err, h_err;
	double			side, A, v_p, mu_0, m, k, tau, g, L, rho, a, b;
	
	side = 1;					// meters
	//side = 1;					// meters
	A = pow(side, 2);				// meters^2, aka 10km^2
	v_p = 5.0/(100*365.25*86400);	// meters/sec, aka 5 cm/year
	mu_0 = 0.05;					// sliding coefficient of friction
	rho = 2.5e3;					// kg/m^3
	m = rho*pow(side, 3);			// kg
	k = sim.G()*sqrt(A)/2;				// spring constant, N/m
	tau = 2*M_PI*sqrt(m/k);			// natural period of vibration
	g = 9.8;						// meters/sec^2
	L = 1e-12;							// meters
	a = 0.003125;
	b = 0.00625;
	
	param_a = a/mu_0;
	param_b = b/mu_0;
	param_k = m*g*mu_0/(k*v_p);
	param_r = pow(tau*v_p/(2*M_PI*L), 2);
	
	//param_a = 0.0625;
	//param_b = 0.125;
	//param_k = 20;
	//param_r = 1e-5;
	std::cerr << param_a << " " << param_b << " " << param_k << " " << param_r << std::endl;
	
	for (i=0;i<NBLOCKS;++i) {
        x_err = v_err = h_err = RCONST(1e-6);
        x_0 = -14.5 + i;
        h_0 = 1;
        v_0 = 1;
		BlockData	bdata(param_a, param_b, param_k, param_r, x_0, v_0, h_0, x_err, v_err, h_err);
		sim.add_block(i, bdata);
	}
	
	// Whether to use simple (Dieterich-style) equations or full rate/state equations
	SimEquations		*eqns;
	
#ifdef USE_SIMPLIFIED_EQNS
	eqns = new SimpleEqns;
#else
	eqns = new OrigEqns;
#endif
	
	CVODESolver *solver = new CVODESolver(eqns);
	//RK4Solver *solver = new RK4Solver(eqns);
	
	// Set the threshold for a rupture to be 0.1 m/s
	sim.set_rupture_threshold(0.1);
	
	// Set the timesteps for each solver (in seconds)
	solver->set_timesteps(1, 0.1);
	solver->set_timesteps(86400, 0.5);
	
	res = sim.init(solver);
	if (res) {
		std::cerr << "Initializing error." << std::endl;
		exit(-1);
	}
	
	fp = fopen("out.txt", "w");
	sim.write_header(fp);
	while(sim.get_time() <= 1000) {
		res = sim.advance();
		if (res != 0) {
			std::cerr << "Err " << res << std::endl;
			sim.write_summary_header(stderr);
			sim.write_summary(stderr);
			break;
		}
		sim.write_cur_data(fp);
	}
	
	fclose(fp);
	
	sim.write_summary_header(stderr);
	sim.write_summary(stderr);
	
	sim.solver()->print_stats();
	
	sim.cleanup();
}
