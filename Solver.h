#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */
#include <cvode/cvode_spgmr.h>
#include <cvode/cvode_spbcgs.h>

#include "Equations.h"

#ifndef _VICARS_SOLVER_H_
#define _VICARS_SOLVER_H_

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

class Solver {
protected:
	SimEquations				*_eqns;
    double                      _timestep;
    
public:
    Solver(SimEquations *eqns) : _eqns(eqns) {};
    virtual ~Solver(void) {};
	virtual int init_solver(ViCaRS *sim) = 0;
	virtual int reinit_solver(realtype cur_time, N_Vector vars) = 0;
	virtual int advance(ViCaRS *sim, N_Vector vars, realtype target_time, realtype &finish_time, bool &next_solver) = 0;
	virtual void print_stats(void) = 0;
    
    void set_timestep(double new_timestep) { _timestep = new_timestep; };
};

class CVODESolver : public Solver {
private:
	// Pointer to the relevant solver
	void					*_solver;
	
	// Vector of initial values, absolute tolerances, and calculated values arranged as follows:
	// [Block1 X, Block1 V, Block1 H, Block2 X, Block2 X, Block2 H, ...]
	N_Vector				_abs_tol;
    
	realtype				_rel_tol;
	
	// Statistics on number of solver steps
	SolverStats             _stats;
    
    int                     _rootdir;
    
    int init_cvode_solver(void **created_solver, int rootdir, ViCaRS *sim);
    
public:
    CVODESolver(SimEquations *eqns) : Solver(eqns), _solver(NULL), _rootdir(0) {};
    ~CVODESolver(void);
	virtual int reinit_solver(realtype cur_time, N_Vector vars) {
        int flag;
        flag = CVodeReInit(_solver, cur_time, vars);
        if (flag != CV_SUCCESS) return 1;
        return 0;
    };
    virtual void set_rootdir(int new_rootdir) { _rootdir = new_rootdir; };
	virtual int init_solver(ViCaRS *sim);
	virtual int advance(ViCaRS *sim, N_Vector vars, realtype target_time, realtype &finish_time, bool &next_solver);
    
	void update_stats(void *solver, SolverStats &stats);
	virtual void print_stats(void);
};

class RK4Solver : public Solver {
private:
    unsigned int        _num_evals;
    N_Vector            k1, k2, k3, k4, in_vals, subsum;
    
public:
    RK4Solver(SimEquations *eqns) : Solver(eqns) {};
	virtual int init_solver(ViCaRS *sim);
	virtual int reinit_solver(realtype cur_time, N_Vector vars) {};
	virtual int advance(ViCaRS *sim, N_Vector vars, realtype target_time, realtype &finish_time, bool &next_solver);
	virtual void print_stats(void);
};

#endif
