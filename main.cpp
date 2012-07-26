#include "RateState.h"

#define NBLOCKS			9

int main(int argc, char **argv)
{
	ViCaRS			sim(NBLOCKS);
	unsigned int	i;
	double			param_a, param_b, param_k, param_r;
	FILE			*fp;
	int				res;
	
	MPI_Init(&argc, &argv);
	
	param_a = 0.0625;
	param_b = 0.125;
	param_k = 20;
	param_r = 1e-5;
	
	for (i=0;i<NBLOCKS;++i) {
		BlockData	bdata(i, param_a, param_b, param_k, param_r, -14.0+i, 1, 1, 1e-8, 1e-8, 1e-8);
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
		if (res != 0) std::cerr << "Err " << res << " t: " << sim.get_time() << std::endl;
		sim.write_cur_data(fp);
	}

	fclose(fp);
	
	sim.print_stats();
	
	sim.cleanup();
	MPI_Finalize();
}
