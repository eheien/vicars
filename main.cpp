
#include "RateState.h"

#define NBLOCKS			1

int main(int argc, char **argv)
{
	ViCaRS			sim(NBLOCKS);
	unsigned int	i;
	double			param_a, param_b, param_k, param_r;
	FILE			*fp;
	int				res;
	
	realtype		x_0, v_0, h_0, x_err, v_err, h_err;
	realtype		sim_end;
	double			A, v_p, mu_0, m, k, tau, g, rho;
	
    // Set global simulation parameters
    sim.params().V_star = 1.0/(1.0e2*1.0e3*365.25*86400);	// 1 centimeter/1000 years in m/s
    sim.params().G = 1.0e4*1.0e6;            // Pascals
    sim.params().D_c = 0.01*1e-3;       // meters
    sim.params().W = 6000;              // meters
	sim.params().beta = 5000;           // meters/second
    sim.params().sigma = 15e6;          // Pascals
    sim.params().A = 0.003;
    sim.params().B = 0.015;
    sim.params().v_ss = 1.0/(1.0e2*365.25*86400);   // 1 cm/year in meters/second
    sim.params().mu_0 = 0.5;            // nominal coefficient of friction
    sim.params().side = 1e2;            // 1 m on a side
    sim.params().L = 1e-12;

	A = pow(sim.params().side, 2);				// meters^2, aka 10km^2
	v_p = sim.params().v_ss;        // meters/sec
	mu_0 = sim.params().mu_0;		// sliding coefficient of friction
	rho = 2.5e3;					// kg/m^3
	m = rho*pow(sim.params().side, 3);			// kg
	k = sim.params().G*sim.params().side/2;	// spring constant, N/m
	tau = 2*M_PI*sqrt(m/k);			// natural period of vibration
	g = 9.8;						// meters/sec^2
	
	param_a = sim.params().A/mu_0;
	param_b = sim.params().B/mu_0;
	param_k = m*g*mu_0/(k*v_p);         // kg*(m/s^2)*(N/m^2)/((N/m)*(m/s))=kg*m/s
	param_r = pow(tau*v_p/(2*M_PI*sim.params().L), 2);   // (sqrt(kg/(N/m))*(m/s)/m)^2=unitless
	
	std::cerr << param_a << " " << param_b << " " << param_k << " " << param_r << std::endl;
	param_a = 0.0625;
	param_b = 0.125;
	//param_k = 20;
	//param_r = 1e-2;
	std::cerr << param_a << " " << param_b << " " << param_k << " " << param_r << std::endl;
    
	for (i=0;i<NBLOCKS;++i) {
        x_err = RCONST(1e-2);
        v_err = RCONST(1e-10);
        h_err = RCONST(1e-5);
        x_0 = -14.5 + i;
        h_0 = 1;
        v_0 = 1;
		BlockData	bdata(param_a, param_b, param_k, param_r, x_0, v_0, h_0, x_err, v_err, h_err);
		sim.add_block(i, bdata);
	}
	
	// Whether to use simple (Dieterich-style) equations or full rate/state equations
	SimEquations		*eqns;
	bool				use_simplified = false, use_cvode = true;
	
	if (use_simplified) {
		eqns = new SimpleEqns;
		sim_end = 50*365.25*86400;
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
		sim_end = 1000*86400;      // 1000 normalized days
		if (use_cvode) {
			CVODESolver *solver_long = new CVODESolver(eqns);
			CVODESolver *solver_rupture = new CVODESolver(eqns);
            realtype long_timestep = 3600;    // 1 day normalized long timesteps
            realtype rupture_timestep = 1;     // 1 second normalized rupture timesteps
            std::cerr << "End: " << sim_end << " Long DT: " << long_timestep << " Short DT: " << rupture_timestep << std::endl;
			solver_long->set_timestep(long_timestep);
			solver_rupture->set_timestep(rupture_timestep);
            solver_long->set_rootdir(1);
            solver_rupture->set_rootdir(-1);
			
			sim.add_solver(solver_rupture);
			sim.add_solver(solver_long);
		} else {
			// RK4
		}
	}
	
	// Set the threshold for a rupture to be v_ss m/s
	sim.set_rupture_threshold(0.001);
	
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
