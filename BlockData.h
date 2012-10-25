#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */

#include <map>

#ifndef _VICARS_BLOCKDATA_H_
#define _VICARS_BLOCKDATA_H_

typedef unsigned int BlockGID;

class BlockData {
public:
	realtype	_a, _b, _k, _r;
	
	BlockData(void) : _a(0), _b(0), _k(0), _r(0) {};
};

typedef std::map<BlockGID, BlockData> BlockMap;


#endif
