
#include "RateState.h"

#define NBLOCKS			7

int main(int argc, char **argv)
{
	ViCaRS			sim(NBLOCKS);
	unsigned int	i;
	FILE			*fp;
	int				res;
	
	realtype		sim_end;
	
    // Set global simulation parameters
    sim.params().V_star = 1.0/(1.0e2*1.0e3*365.25*86400);	// 1 centimeter/1000 years in m/s
    sim.params().G = 3.0e10;            // Pascals
    sim.params().D_c = 0.01*1e-3;       // meters
    sim.params().W = 6000;              // meters
	sim.params().beta = 5000;           // meters/second
    sim.params().sigma = 15e6;          // Pascals
    //sim.params().A = 0.003;
    //sim.params().B = 0.015;
    sim.params().A = 0.003125;
    sim.params().B = 0.00625;
    sim.params().v_ss = 10.0/(1e2*365.25*86400);   // 1 m/year in meters/second
    sim.params().mu_0 = 0.05;            // nominal coefficient of friction
    //sim.params().mu_0 = 0.5;            // nominal coefficient of friction
    sim.params().side = 8;            // 1 km on a side
    sim.params().rho = 2.5e3;           // kg/m^3 of rock
    sim.params().gravity = 9.8;         // m/s^2
    sim.params().L = 1e-6;

	for (i=0;i<NBLOCKS;++i) {
		BlockData	bdata;
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
        realtype time_norm_factor = sim.params().v_ss/sim.params().L;
        //std::cerr << "Norm factor " << time_norm_factor << std::endl;
		sim_end = 1000;      // 1000 normalized days
		if (use_cvode) {
			CVODESolver *solver_long = new CVODESolver(eqns);
			CVODESolver *solver_rupture = new CVODESolver(eqns);
            realtype long_timestep = 0.05;    // 1 day normalized long timesteps
            realtype rupture_timestep = 0.01;     // 1 second normalized rupture timesteps
            //std::cerr << "End: " << sim_end << " Long DT: " << long_timestep << " Short DT: " << rupture_timestep << std::endl;
			solver_long->set_timestep(long_timestep);
			solver_rupture->set_timestep(rupture_timestep);
            solver_long->set_rupture_threshold(0.01);
            solver_rupture->set_rupture_threshold(0.001);
            solver_long->set_rootdir(1);
            solver_rupture->set_rootdir(-1);
			
			sim.add_solver(solver_long);
			//sim.add_solver(solver_rupture);
		} else {
			// RK4
		}
	}
	
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
        fflush(fp);
	}
	
	fclose(fp);
	
	sim.write_summary_header(stderr);
	sim.write_summary(stderr);
	
	//sim.solver()->print_stats();
	
	sim.cleanup();
}
