#include "RateState.h"

#define NBLOCKS			1

int main(int argc, char **argv)
{
	ViCaRS			sim(NBLOCKS);
	unsigned int	i;
	double			param_a, param_b, param_k, param_r;
	FILE			*fp;
	int				res, world_size, rank;
	
	realtype		x_0, v_0, h_0, x_err, v_err, h_err;
	double			side, A, G, v_p, mu_0, m, k, tau, g, L, rho, a, b;
	
	MPI_Init(&argc, &argv);
	
	side = 10*1000;					// meters
	//side = 1;					// meters
	A = pow(side, 2);				// meters^2, aka 10km^2
	G = 3.0e10;						// Pascals
	v_p = 5.0/(100*365.25*86400);	// meters/sec, aka 5 cm/year
	mu_0 = 0.05;					// sliding coefficient of friction
	rho = 2.5e3;					// kg/m^3
	m = rho*pow(side, 3);			// kg
	k = G*sqrt(A)/2;				// spring constant, N/m
	tau = 2*M_PI*sqrt(m/k);			// natural period of vibration
	g = 9.8;						// meters/sec^2
	L = 1;							// meters
	a = 0.003125;
	b = 0.00625;
	
	param_a = a/mu_0;
	param_b = b/mu_0;
	param_k = m*g*mu_0/(k*v_p);
	param_r = pow(tau*v_p/(2*M_PI*L), 2);
	
	param_a = 0.0625;
	param_b = 0.125;
	param_k = 20;
	param_r = 1e-5;
	
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	for (i=0;i<NBLOCKS;++i) {
        x_err = v_err = h_err = RCONST(1e-6);
        x_0 = -14.5 + i;
        h_0 = 1;
        v_0 = 1;
		BlockData	bdata(i, param_a, param_b, param_k, param_r, x_0, v_0, h_0, x_err, v_err, h_err);
		sim.add_local_block(bdata);
	}
	
	// Set the threshold for a rupture to be 0.1 m/s
	sim.set_rupture_threshold(0.1);
	
	// Set the timesteps for each solver (in seconds)
	sim.set_timesteps(1, 0.1);
	
	res = sim.init();
	if (res) {
		std::cerr << "Initializing error." << std::endl;
		exit(-1);
	}
	
	fp = fopen("out.txt", "w");
	sim.write_header(fp);
	while(sim.get_time() <= 500) {
		res = sim.advance();
		if (res != 0) {
			std::cerr << "Err " << res << std::endl;
			sim.write_summary_header(stderr);
			sim.write_summary(stderr);
			break;
		}
		//else std::cout << "successfully advanced (t: " << sim.get_time() << ")" << std::endl;
		sim.write_cur_data(fp);
	}
	
	fclose(fp);
	
	sim.print_stats();
	
	sim.cleanup();
	MPI_Finalize();
}
