//#include "QuakeLib.h"
#include "Spline.h"
#include "Equations.h"
#include "Solver.h"
#include "Params.h"

#include <stdio.h>
#include <math.h>
#include <vector>
#include <map>
#include <iostream>

#ifndef _VICARS_RATESTATE_H_
#define _VICARS_RATESTATE_H_

// Format of result vector: [x0, v0, h0, x1, v1, h1, ..., xn, vn, hn]

class ViCaRS {
private:
	unsigned int			_num_global_blocks;
	unsigned int			_num_equations;
	
    SimParams               _params;
    
	// Vector of initial values, absolute tolerances, and calculated values arranged as follows:
	// [Block1 X, Block1 V, Block1 H, Block2 X, Block2 X, Block2 H, ...]
	N_Vector				_vars;
	
	BlockMap				_bdata;
	
	realtype				_cur_time;
	
	// Generic solver class
	std::vector<Solver*>	_solvers;
	unsigned int			_cur_solver;
	
	SimEquations			*_eqns;
	
	// Whether to use the slowness law (var = true) or slip law (var = false)
	bool                    _use_slowness_law;
	
	//quakelib::DenseStdStraight<double> greens_matrix;
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
	
	N_Vector vars(void) const { return _vars; };
	
	int add_block(const BlockGID &id, const BlockData &block_data);
	
    SimParams &params(void) { return _params; };
    
	int init(SimEquations *eqns);
	int advance(void);
	void cleanup(void);
	void add_solver(Solver *solver) { _solvers.push_back(solver); };
	SimEquations *equations(void) { return _eqns; };
	
	realtype interaction(BlockGID a, BlockGID b);
	realtype get_time(void) const { return _cur_time; };
	
    int cur_solver(void) const { return _cur_solver; };
    realtype rupture_threshold(void) { return _solvers[_cur_solver]->rupture_threshold(); };
    
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
};

#endif
