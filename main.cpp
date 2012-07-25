#include "RateState.h"

#define NBLOCKS			1

int main(int argc, char **argv)
{
	ViCaRS			sim(NBLOCKS);
	unsigned int	i;
	double			param_a, param_b, param_k, param_r;
	realtype		t=0, tstep=0.5;
	FILE			*fp;
	int				res;
	
	MPI_Init(&argc, &argv);
	
	param_a = 0.0625;
	param_b = 0.125;
	param_k = 20;
	param_r = 1e-5;
	
	for (i=0;i<NBLOCKS;++i) {
		BlockData	bdata(i, param_a, param_b, param_k, param_r, -14.0+i, 1, 1, 1e-6, 1e-6, 1e-6);
		sim.add_local_block(bdata);
	}
	
	res = sim.init();
	if (res) {
		std::cerr << "Initializing error." << std::endl;
		exit(-1);
	}
	
	fp = fopen("out.txt", "w");
	sim.write_header(fp);
	t = tstep;
	while(t <= 500) {
		res = sim.advance(t);
		if (res != 0) std::cerr << "Err " << res << " t: " << t << std::endl;
		sim.write_cur_data(fp);
		t += tstep;
	}

	fclose(fp);
	
	sim.print_stats();
	
	sim.cleanup();
	MPI_Finalize();
}
