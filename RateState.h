#include "Matrix.h"
#include "Spline.h"
#include "Equations.h"
#include "Solver.h"

#include <stdio.h>
#include <math.h>
#include <vector>
#include <map>
#include <iostream>

#ifndef _VICARS_RATESTATE_H_
#define _VICARS_RATESTATE_H_

// Format of result vector: [x0, v0, h0, x1, v1, h1, ..., xn, vn, hn]

/* Problem Constants */

int solve_odes(realtype t, N_Vector y, N_Vector ydot, void *user_data);
int jacobian_times_vector(N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, void *user_data, N_Vector tmp);
int check_for_mode_change(realtype t, N_Vector y, realtype *gout, void *user_data);

class SolverStats {
public:
	long int            _nst, _nfe, _nni, _ncfn, _netf, _njtv, _nge;
    
    SolverStats(void) : _nst(0), _nfe(0), _nni(0), _ncfn(0), _netf(0), _njtv(0), _nge(0) {};
    void add_counts(long int nst, long int nfe, long int nni, long int ncfn, long int netf, long int njtv, long int nge) {
        _nst += nst;
        _nfe += nfe;
        _nni += nni;
        _ncfn += ncfn;
        _netf += netf;
        _njtv += njtv;
        _nge += nge;
    }
};

class ViCaRS {
private:
	unsigned int			_num_global_blocks;
	unsigned int			_num_equations;
	
	// Vector of initial values, absolute tolerances, and calculated values arranged as follows:
	// [Block1 X, Block1 V, Block1 H, Block2 X, Block2 X, Block2 H, ...]
	N_Vector				_abs_tol;
	N_Vector				_vars;
	
	BlockMap				_bdata;
	
	realtype				_G;		// shear modulus parameter
	
	realtype				_rel_tol;
	realtype				_cur_time;
	
	// CVODE solvers
	void					*_solver_long, *_solver_rupture;
	
	// Pointer to the currently active solver
	void					*_current_solver;
	
	realtype                _rupture_timestep, _long_timestep;
	
	// Whether the simulation is in rupture mode (var = true) or long term mode (var = false)
	bool                    _in_rupture;
	
	// When any block in the system reaches this speed, the system is considered to be in rupture mode
	realtype                _rupture_threshold;
	
	// Statistics on number of solver steps for long term and rupture solvers
	SolverStats             _stats_long, _stats_rupture;
	
	// Whether to use the slowness law (var = true) or slip law (var = false)
	bool                    _use_slowness_law;
	
	EqnSolver				*_eqns;
	
	VCDenseStdStraight greens_matrix;
	int fill_greens_matrix(void);
	
public:
	typedef BlockMap::iterator iterator;
	typedef BlockMap::const_iterator const_iterator;
	
	ViCaRS(unsigned int total_num_blocks);
	~ViCaRS(void) {};
	
	const_iterator begin(void) const { return _bdata.begin(); };
	iterator begin(void) { return _bdata.begin(); };
	const_iterator end(void) const { return _bdata.end(); };
	iterator end(void) { return _bdata.end(); };
	
	int add_block(const BlockGID &id, const BlockData &block_data);
	
	int init(EqnSolver *eqns);
	int init_solver(void **solver, int rootdir);
	int advance(void);
	int advance_simple(void);
	void cleanup(void);
	
	EqnSolver *get_eqns(void) { return _eqns; };
	
	realtype G(void) const { return _G; };
	
	realtype interaction(BlockGID a, BlockGID b);
	void set_timesteps(realtype long_term_step, realtype rupture_step) { _long_timestep = long_term_step; _rupture_timestep = rupture_step; };
	realtype get_time(void) const { return _cur_time; };
	
	void set_rupture_threshold(realtype new_threshold) { _rupture_threshold = new_threshold; };
	realtype rupture_threshold(void) const { return _rupture_threshold; };
	
	bool use_slowness_law(void) const { return _use_slowness_law; };
	
	unsigned int num_global_blocks(void) const { return _num_global_blocks; };
	unsigned int num_eqs(void) const { return _num_equations; };
	
	realtype &param_a(BlockGID block_num) { return _bdata[block_num]._a; };
	realtype &param_b(BlockGID block_num) { return _bdata[block_num]._b; };
	realtype &param_k(BlockGID block_num) { return _bdata[block_num]._k; };
	realtype &param_r(BlockGID block_num) { return _bdata[block_num]._r; };
	
	void write_header(FILE *fp);
	void write_summary(FILE *fp);
	void write_summary_header(FILE *fp);
	void write_cur_data(FILE *fp);
	
	void print_stats(void);
	void update_stats(void *solver, SolverStats &stats);
};

#endif
