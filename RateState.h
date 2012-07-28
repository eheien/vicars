#include <stdio.h>
#include <math.h>
#include <vector>
#include <iostream>
#include "Matrix.h"
#include "Spline.h"

#include <mpi.h>

#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>
#include <nvector/nvector_parallel.h>
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */
#include <cvode/cvode_spgmr.h>
#include <cvode/cvode_spbcgs.h>

// Format of result vector: [x0, v0, h0, x1, v1, h1, ..., xn, vn, hn]

#define Xth(y,block)	NV_Ith_P(y,block*NEQ+0)
#define Vth(y,block)	NV_Ith_P(y,block*NEQ+1)
#define Hth(y,block)	NV_Ith_P(y,block*NEQ+2)

/* Problem Constants */

#define NEQ				3
#define EQ_X			0
#define EQ_V			1
#define EQ_H			2

int func(realtype t, N_Vector y, N_Vector ydot, void *user_data);
int rupture_test(realtype t, N_Vector y, realtype *gout, void *user_data);
int jacobian_times_vector(N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, void *user_data, N_Vector tmp);

int simple_equations(realtype t, N_Vector y, N_Vector ydot, void *user_data);
int simple_rupture_test(realtype t, N_Vector y, realtype *gout, void *user_data);
int simple_jacobian_times_vector(N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, void *user_data, N_Vector tmp);

typedef unsigned int BlockGID;

typedef std::map<BlockGID, unsigned int> GlobalLocalMap;

class BlockData {
public:
	BlockGID	_gid;
	realtype	_a, _b, _k, _r;
	realtype	_init_x, _init_v, _init_h;
	realtype	_tol_x, _tol_v, _tol_h;
	
	BlockData(BlockGID gid, realtype a, realtype b, realtype k, realtype r,
			  realtype init_x, realtype init_v, realtype init_h,
			  realtype tol_x, realtype tol_v, realtype tol_h) :
				_gid(gid), _a(a), _b(b), _k(k), _r(r),
				_init_x(init_x), _init_v(init_v), _init_h(init_h),
				_tol_x(tol_x), _tol_v(tol_v), _tol_h(tol_h) {};
};

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
	GlobalLocalMap		_global_local_map;

	// Vector of initial values, absolute tolerances, and calculated values arranged as follows:
	// [Block1 X, Block1 V, Block1 H, Block2 X, Block2 X, Block2 H, ...]
	N_Vector				_abs_tol;
	N_Vector				_vars;
	std::vector<BlockData>	_bdata;
	
	realtype				_rel_tol;
	realtype				_cur_time;

	// CVODE solvers
	void					*_solver_long, *_solver_rupture;
	// Pointer to the currently active solver
	void					*_current_solver;
	
	// Communication related variables
	int						world_size, rank;
	
	realtype                _rupture_timestep, _long_timestep;
    
	// Whether the simulation is in rupture mode (var = true) or long term mode (var = false)
  bool                    _in_rupture;
    
	// When any block in the system reaches this speed, the system is considered to be in rupture mode
  realtype                _rupture_threshold;
    
  // Statistics on number of solver steps for long term and rupture solvers
  SolverStats             _stats_long, _stats_rupture;
    
  // Whether to use the slowness law (var = true) or slip law (var = false)
  bool                    _use_slowness_law;
    
	LogSpline	            log_approx;
    
  // Whether to use the log spline approximation or not
  bool                    _use_log_spline;
   
  VCDenseStdStraight greens_matrix;
  int fill_greens_matrix(void);
 
public:
	typedef GlobalLocalMap::iterator iterator;
	typedef GlobalLocalMap::const_iterator const_iterator;
	
	ViCaRS(unsigned int total_num_blocks);
	~ViCaRS(void) {};
	
	const_iterator begin(void) const { return _global_local_map.begin(); };
	iterator begin(void) { return _global_local_map.begin(); };
	const_iterator end(void) const { return _global_local_map.end(); };
	iterator end(void) { return _global_local_map.end(); };
	
	int add_local_block(const BlockData &block_data);
	unsigned int global_to_local(const BlockGID &gid) { return _global_local_map[gid]; };
	
	int init(void);
	int advance(void);
	void cleanup(void);
	
	realtype interaction(BlockGID a, BlockGID b);
	void set_timesteps(realtype long_term_step, realtype rupture_step) { _long_timestep = long_term_step; _rupture_timestep = rupture_step; };
	realtype get_time(void) const { return _cur_time; };
    
	void set_rupture_threshold(realtype new_threshold) { _rupture_threshold = new_threshold; };
	realtype rupture_threshold(void) const { return _rupture_threshold; };
	
    bool use_slowness_law(void) const { return _use_slowness_law; };
    bool use_log_spline(void) const { return _use_log_spline; };
    
	unsigned int num_global_blocks(void) const { return _num_global_blocks; };
	unsigned int num_eqs(void) const { return _num_equations; };
	
	realtype &param_a(BlockGID block_num) { unsigned int lid=_global_local_map[block_num]; return _bdata[lid]._a; };
	realtype &param_b(BlockGID block_num) { unsigned int lid=_global_local_map[block_num]; return _bdata[lid]._b; };
	realtype &param_k(BlockGID block_num) { unsigned int lid=_global_local_map[block_num]; return _bdata[lid]._k; };
	realtype &param_r(BlockGID block_num) { unsigned int lid=_global_local_map[block_num]; return _bdata[lid]._r; };
	
	realtype &X(BlockGID block_num) { unsigned int lid=_global_local_map[block_num]; return NV_DATA_P(_vars)[lid*_num_equations+0]; };
	realtype &V(BlockGID block_num) { unsigned int lid=_global_local_map[block_num]; return NV_DATA_P(_vars)[lid*_num_equations+1]; };
	realtype &H(BlockGID block_num) { unsigned int lid=_global_local_map[block_num]; return NV_DATA_P(_vars)[lid*_num_equations+2]; };
	
    realtype log_func(realtype x) const {
        if (_use_log_spline) return log_approx(x);
        else return log(x);
    }
	
    realtype log_deriv(realtype x) const {
        if (_use_log_spline) return log_approx.deriv(x);
        else return 1/x;
    }
    
	realtype F(BlockGID block_num, realtype v, realtype h) {
		return 1 + param_a(block_num) * log_func(v) + param_b(block_num)*log_func(h);
	}
	
	void write_header(FILE *fp);
	void write_cur_data(FILE *fp);
	
	void print_stats(void);
	void update_stats(void *solver, SolverStats &stats);
};
