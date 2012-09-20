#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */

#include <map>

#ifndef _VICARS_BLOCKDATA_H_
#define _VICARS_BLOCKDATA_H_

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


#endif
