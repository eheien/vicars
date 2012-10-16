
#include "RateState.h"

#define NBLOCKS			1
#define USE_CVODE
#define USE_SIMPLIFIED_EQNS

int main(int argc, char **argv)
{
	ViCaRS			sim(NBLOCKS);
	unsigned int	i;
	double			param_a, param_b, param_k, param_r;
	FILE			*fp;
	int				res;
	
	realtype		x_0, v_0, h_0, x_err, v_err, h_err;
	realtype		sim_end;
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
	param_k = 20;
	param_r = 1e-5;
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
	bool				use_simplified = true, use_cvode = true;
	
	if (use_simplified) {
		eqns = new SimpleEqns;
		sim_end = 10*365.25*86400;
		if (use_cvode) {
			CVODESolver *solver_phase0 = new CVODESolver(eqns);
			CVODESolver *solver_phase1 = new CVODESolver(eqns);
			CVODESolver *solver_phase2 = new CVODESolver(eqns);
			solver_phase0->set_timestep(86400);
			solver_phase1->set_timestep(60);
			solver_phase2->set_timestep(1);
			
			sim.add_solver(solver_phase0);
			sim.add_solver(solver_phase1);
			sim.add_solver(solver_phase2);
		} else {
			// RK4
		}
	} else {
		eqns = new OrigEqns;
		sim_end = 1000;
		if (use_cvode) {
			CVODESolver *solver_long = new CVODESolver(eqns);
			CVODESolver *solver_rupture = new CVODESolver(eqns);
			solver_long->set_timestep(1);
			solver_rupture->set_timestep(0.1);
			
			sim.add_solver(solver_long);
			sim.add_solver(solver_rupture);
		} else {
			// RK4
		}
	}
	
	// Set the threshold for a rupture to be 0.1 m/s
	sim.set_rupture_threshold(0.1);
	
	res = sim.init(eqns);
	if (res) {
		std::cerr << "Initializing error." << std::endl;
		exit(-1);
	}
	
	fp = fopen("out.txt", "w");
	sim.write_header(fp);
	while(sim.get_time() <= sim_end) {
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
	
	//sim.solver()->print_stats();
	
	sim.cleanup();
}
