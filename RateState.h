#include <stdio.h>
#include <math.h>
#include <vector>
#include <map>
#include <iostream>
#include "Matrix.h"
#include "Spline.h"

#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */
#include <cvode/cvode_spgmr.h>
#include <cvode/cvode_spbcgs.h>

// Format of result vector: [x0, v0, h0, x1, v1, h1, ..., xn, vn, hn]

#define Xth(y,block)	NV_Ith_S(y,block*NEQ+0)
#define Vth(y,block)	NV_Ith_S(y,block*NEQ+1)
#define Hth(y,block)	NV_Ith_S(y,block*NEQ+2)

/* Problem Constants */

#define NEQ				3
#define EQ_X			0
#define EQ_V			1
#define EQ_H			2

int solve_odes(realtype t, N_Vector y, N_Vector ydot, void *user_data);
int jacobian_times_vector(N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, void *user_data, N_Vector tmp);
int check_for_rupture(realtype t, N_Vector y, realtype *gout, void *user_data);

typedef unsigned int BlockGID;

class BlockData {
public:
	realtype	_a, _b, _k, _r;
	realtype	_init_x, _init_v, _init_h;
	realtype	_tol_x, _tol_v, _tol_h;
	
	BlockData(void) : _a(0), _b(0), _k(0), _r(0), _init_x(0), _init_v(0), _init_h(0), _tol_x(0), _tol_v(0), _tol_h(0) {};
	BlockData(realtype a, realtype b, realtype k, realtype r,
			  realtype init_x, realtype init_v, realtype init_h,
			  realtype tol_x, realtype tol_v, realtype tol_h) :
	_a(a), _b(b), _k(k), _r(r),
	_init_x(init_x), _init_v(init_v), _init_h(init_h),
	_tol_x(tol_x), _tol_v(tol_v), _tol_h(tol_h) {};
};

typedef std::map<BlockGID, BlockData> BlockMap;

class ViCaRS;

class EqnSolver {
public:
	virtual ~EqnSolver(void) {};
	virtual int init(ViCaRS *sim) = 0;
	
	virtual int solve_odes(ViCaRS *sim, realtype t, N_Vector y, N_Vector ydot) = 0;
	virtual bool has_jacobian(void) = 0;
	virtual int jacobian_times_vector(ViCaRS *sim, N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, N_Vector tmp) = 0;
	virtual int check_for_rupture(ViCaRS *sim, realtype t, N_Vector y, realtype *gout) = 0;
	virtual bool values_valid(ViCaRS *sim, N_Vector y);
};

class OrigEqns : public EqnSolver {
public:
	virtual ~OrigEqns(void) {};
	virtual int init(ViCaRS *sim);
	
	virtual int solve_odes(ViCaRS *sim, realtype t, N_Vector y, N_Vector ydot);
	virtual bool has_jacobian(void) { return true; };
	virtual int jacobian_times_vector(ViCaRS *sim, N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, N_Vector tmp);
	virtual int check_for_rupture(ViCaRS *sim, realtype t, N_Vector y, realtype *gout);
};

class SimpleEqns : public EqnSolver {
private:
	// Phase of each element
	std::map<BlockGID, int> phase;
	
	realtype		mu_0, A, B;
	
	N_Vector		_ss_stress, _stress_loading;
	
public:
	SimpleEqns(void) : mu_0(0.5), A(0.005), B(0.015), _ss_stress(NULL), _stress_loading(NULL) {};
	virtual ~SimpleEqns(void) { if (_ss_stress) N_VDestroy_Serial(_ss_stress); if (_stress_loading) N_VDestroy_Serial(_stress_loading); };
	virtual int init(ViCaRS *sim);
	
	virtual int solve_odes(ViCaRS *sim, realtype t, N_Vector y, N_Vector ydot);
	virtual bool has_jacobian(void) { return false; };
	virtual int jacobian_times_vector(ViCaRS *sim, N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, N_Vector tmp) { return -1; };
	virtual int check_for_rupture(ViCaRS *sim, realtype t, N_Vector y, realtype *gout);
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
	
	// Whether to use simple (dieterich-style) equations or full rate/state equations
	bool                    _use_simple_equations;
	
	LogSpline				log_approx;
	
	// Whether to use the log spline approximation or not
	bool                    _use_log_spline;
	
	EqnSolver				*_eqns;
	
	VCDenseStdStraight greens_matrix;
	int fill_greens_matrix(void);
	
public:
	N_Vector _stress;
	realtype h_ss, v_ss, v_eq, v_min;
	
	typedef BlockMap::iterator iterator;
	typedef BlockMap::const_iterator const_iterator;
	
	ViCaRS(unsigned int total_num_blocks);
	~ViCaRS(void) {};
	
	const_iterator begin(void) const { return _bdata.begin(); };
	iterator begin(void) { return _bdata.begin(); };
	const_iterator end(void) const { return _bdata.end(); };
	iterator end(void) { return _bdata.end(); };
	
	int add_block(const BlockGID &id, const BlockData &block_data);
	
	int init(void);
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
	bool use_log_spline(void) const { return _use_log_spline; };
	bool use_simple_equations(void) const { return _use_simple_equations; };
	
	unsigned int num_global_blocks(void) const { return _num_global_blocks; };
	unsigned int num_eqs(void) const { return _num_equations; };
	
	realtype &param_a(BlockGID block_num) { return _bdata[block_num]._a; };
	realtype &param_b(BlockGID block_num) { return _bdata[block_num]._b; };
	realtype &param_k(BlockGID block_num) { return _bdata[block_num]._k; };
	realtype &param_r(BlockGID block_num) { return _bdata[block_num]._r; };
	
	realtype &X(BlockGID block_num) { return NV_DATA_S(_vars)[block_num*_num_equations+0]; };
	realtype &V(BlockGID block_num) { return NV_DATA_S(_vars)[block_num*_num_equations+1]; };
	realtype &H(BlockGID block_num) { return NV_DATA_S(_vars)[block_num*_num_equations+2]; };
	
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
	void write_summary(FILE *fp);
	void write_summary_header(FILE *fp);
	void write_cur_data(FILE *fp);
	
	void print_stats(void);
	void update_stats(void *solver, SolverStats &stats);
};
