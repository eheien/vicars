#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */

#ifndef _VICARS_SOLVER_H_
#define _VICARS_SOLVER_H_

class Solver {
public:
	virtual void init_solver(void) = 0;
	virtual void advance(realtype target_time, realtype actual_time) = 0;
	virtual void cleanup(void) = 0;
};

class CVODESolver : public Solver {
	virtual void init_solver(void);
	virtual void advance(realtype target_time, realtype actual_time);
	virtual void cleanup(void);
};

class RK4Solver : public Solver {
	virtual void init_solver(void);
	virtual void advance(realtype target_time, realtype actual_time);
	virtual void cleanup(void);
};

#endif
