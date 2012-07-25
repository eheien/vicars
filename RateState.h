#include <stdio.h>
#include <math.h>
#include <vector>
#include <iostream>
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

#define USE_LOG_SPLINE

int func(realtype t, N_Vector y, N_Vector ydot, void *user_data);
int jacobian_times_vector(N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, void *user_data, N_Vector tmp);

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

class ViCaRS {
private:
	unsigned int			_num_global_blocks;
	unsigned int			_num_equations;
	GlobalLocalMap			_global_local_map;
	
	// Vector of initial values, absolute tolerances, and calculated values arranged as follows:
	// [Block1 X, Block1 V, Block1 H, Block2 X, Block2 X, Block2 H, ...]
	N_Vector				_abs_tol;
	N_Vector				_vars;
	std::vector<BlockData>	_bdata;
	
	realtype				_rel_tol;
	realtype				_cur_time;

	// CVODE solver
	void					*_solver;
	
	// Communication related variables
	int						world_size, rank;
	int						*num_send, *num_recv;
	int						total_send, total_recv;
	int						*send_offset, *recv_offset;
	
	LogSpline				log_hack;
    
public:
	typedef GlobalLocalMap::iterator iterator;
	typedef GlobalLocalMap::const_iterator const_iterator;
	
	ViCaRS(unsigned int total_num_blocks) :_num_global_blocks(total_num_blocks), _num_equations(NEQ) {};
	
	~ViCaRS(void) {};
	
	const_iterator begin(void) const { return _global_local_map.begin(); };
	iterator begin(void) { return _global_local_map.begin(); };
	const_iterator end(void) const { return _global_local_map.end(); };
	iterator end(void) { return _global_local_map.end(); };
	
	int add_local_block(const BlockData &block_data);
	unsigned int global_to_local(const BlockGID &gid) { return _global_local_map[gid]; };
	
	int init(void);
	int advance(realtype next_time);
	void cleanup(void);
	
	realtype interaction(BlockGID a, BlockGID b);
	
	unsigned int num_global_blocks(void) const { return _num_global_blocks; };
	unsigned int num_eqs(void) const { return _num_equations; };
	
	realtype &param_a(BlockGID block_num) { unsigned int lid=_global_local_map[block_num]; return _bdata[lid]._a; };
	realtype &param_b(BlockGID block_num) { unsigned int lid=_global_local_map[block_num]; return _bdata[lid]._b; };
	realtype &param_k(BlockGID block_num) { unsigned int lid=_global_local_map[block_num]; return _bdata[lid]._k; };
	realtype &param_r(BlockGID block_num) { unsigned int lid=_global_local_map[block_num]; return _bdata[lid]._r; };
	
	//realtype X(BlockGID block_num) const { unsigned int lid=_global_local_map[block_num]; return NV_DATA_P(_vars)[lid*_num_equations+0]; };
	realtype &X(BlockGID block_num) { unsigned int lid=_global_local_map[block_num]; return NV_DATA_P(_vars)[lid*_num_equations+0]; };
	//realtype V(BlockGID block_num) const { unsigned int lid=_global_local_map[block_num]; return NV_DATA_P(_vars)[lid*_num_equations+1]; };
	realtype &V(BlockGID block_num) { unsigned int lid=_global_local_map[block_num]; return NV_DATA_P(_vars)[lid*_num_equations+1]; };
	//realtype H(BlockGID block_num) const { unsigned int lid=_global_local_map[block_num]; return NV_DATA_P(_vars)[lid*_num_equations+2]; };
	realtype &H(BlockGID block_num) { unsigned int lid=_global_local_map[block_num]; return NV_DATA_P(_vars)[lid*_num_equations+2]; };
	
	realtype F(BlockGID block_num, realtype v, realtype h) {
#ifdef USE_LOG_SPLINE
		return 1 + param_a(block_num)*log_hack(v) + param_b(block_num)*log_hack(h);
#else
		return 1 + param_a(block_num)*log(v) + param_b(block_num)*log(h);
#endif
	}
	
	void write_header(FILE *fp);
	void write_cur_data(FILE *fp);
	
	void print_stats(void);
};
