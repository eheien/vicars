#include "Solver.h"
#include "RateState.h"

CVODESolver::~CVODESolver(void) {
	N_VDestroy_Serial(_abs_tol);
    
	// Free integrator memory
	if (_solver_long) CVodeFree(&_solver_long);
	// Free integrator memory
	if (_solver_rupture) CVodeFree(&_solver_rupture);
}

int CVODESolver::init_solver(ViCaRS *sim) {
	int				flag, rootdir;
	BlockMap::const_iterator	it;
	BlockGID		gid;
	
	_in_rupture = false;
	
	// Initialize variables and tolerances
	_rel_tol = RCONST(1.0e-7);
	_abs_tol = N_VNew_Serial(sim->num_global_blocks()*_eqns->num_equations());
	if (_abs_tol == NULL) return 1;
	
	for (it=sim->begin();it!=sim->end();++it) {
		gid = it->first;
		_eqns->init_block(gid, it->second, sim->vars(), _abs_tol);
	}
    
	// Initialize the long term solver
	rootdir = 1;
	flag = init_cvode_solver(&_solver_long, rootdir, sim);
	if (flag) return flag;
	
	// And the rupture solver
	rootdir = -1;
	flag = init_cvode_solver(&_solver_rupture, rootdir, sim);
	if (flag) return flag;
	
	// Set the current solver to be the long term
	_current_solver = _solver_long;
	
	// Init CVodes structures
	flag = _eqns->init(sim);
	if (flag) return flag;
    
    return 0;
}

int CVODESolver::init_cvode_solver(void **created_solver, int rootdir, ViCaRS *sim) {
	void	*solver;
	int		flag;
	
	// Call CVodeCreate to create the solver memory and specify the
	// Backward Differentiation Formula and the use of a Newton iteration
	solver = CVodeCreate(CV_BDF, CV_NEWTON);
	if (solver == NULL) return 1;
	
	// Turn off error messages
	//flag = CVodeSetErrFile(solver, NULL);
	//if (flag != CV_SUCCESS) return 1;
    
	// Call CVodeInit to initialize the integrator memory and specify the
	// user's right hand side function in y'=f(t,y), the inital time, and
	// the initial dependent variable vector.
	flag = CVodeInit(solver, solve_odes, 0, sim->vars());     // TODO: Change 0 to _cur_time
	if (flag != CV_SUCCESS) return 1;
	
	// Call CVodeSVtolerances to specify the scalar relative tolerance
	// and vector absolute tolerances
	flag = CVodeSVtolerances(solver, _rel_tol, _abs_tol);
	if (flag != CV_SUCCESS) return 1;
	
	// Set the root finding function, going above the rupture limit for long term solver
	// and going below the limit for rupture solver.
	flag = CVodeRootInit(solver, 1, check_for_mode_change);
	if (flag != CV_SUCCESS) return 1;
	flag = CVodeSetRootDirection(solver, &rootdir);
	if (flag != CV_SUCCESS) return 1;
	
	// Call CVSpbcg to specify the CVSPBCG scaled preconditioned Bi-CGSTab iterative solver
	// Currently not using any preconditioner
	flag = CVSpbcg(solver, PREC_NONE, 0);
	if (flag != CV_SUCCESS) return 1;
    
	// Set the user data to be the main simulation object
	flag = CVodeSetUserData(solver, sim);
	if (flag != CV_SUCCESS) return 1;
	
	// Set maximum number of steps per solver iteration
	// This should be higher if the time step is less frequent or tolerance is lowered
	flag = CVodeSetMaxNumSteps(solver, 1e7);
	if (flag != CV_SUCCESS) return 1;
	
	// Set the Jacobian x vector function
	if (_eqns->has_jacobian()) {
		flag = CVSpilsSetJacTimesVecFn(solver, jacobian_times_vector);
		if (flag != CV_SUCCESS) return 1;
	}
	
	*created_solver = solver;
	
	return 0;
}

int CVODESolver::advance(ViCaRS *sim, N_Vector vars, realtype target_time, realtype &finish_time) {
	int	flag;
	realtype tstep;
	
	if (_in_rupture) tstep = _rupture_timestep;
	else tstep = _long_timestep;
	
	flag = CVode(_current_solver, target_time+tstep, vars, &finish_time, CV_NORMAL);
	if (flag == CV_SUCCESS) {
		return 0;
	} else if (flag == CV_ROOT_RETURN) {
		// Update the stats for each solver
		update_stats(_solver_long, _stats_long);
		update_stats(_solver_rupture, _stats_rupture);
		
        // Change the current solver
		if (_in_rupture) _current_solver = _solver_long;
		else _current_solver = _solver_rupture;
		
		// Reinit both solvers
		flag = CVodeReInit(_solver_long, finish_time, vars);
		if (flag != CV_SUCCESS) return 1;
		flag = CVodeReInit(_solver_rupture, finish_time, vars);
		if (flag != CV_SUCCESS) return 1;
		_in_rupture = !_in_rupture;
		
		_eqns->handle_mode_change(sim, finish_time, vars);
		
		return 0;
	} else return flag;
}

void CVODESolver::update_stats(void *solver, SolverStats &stats) {
	long int nst, nfe, nni, ncfn, netf, njtv, nge;
	int flag;
	
	flag = CVodeGetNumSteps(solver, &nst);
	if (flag!=CV_SUCCESS) printf("failed to update SolverStats");
	flag = CVodeGetNumRhsEvals(solver, &nfe);
	if (flag!=CV_SUCCESS) printf("failed to update SolverStats");
	flag = CVodeGetNumErrTestFails(solver, &netf);
	if (flag!=CV_SUCCESS) printf("failed to update SolverStats");
	flag = CVodeGetNumNonlinSolvIters(solver, &nni);
	if (flag!=CV_SUCCESS) printf("failed to update SolverStats");
	flag = CVodeGetNumNonlinSolvConvFails(solver, &ncfn);
	if (flag!=CV_SUCCESS) printf("failed to update SolverStats");
	flag = CVSpilsGetNumJtimesEvals(solver, &njtv);
	if (flag!=CV_SUCCESS) printf("failed to update SolverStats");
	flag = CVodeGetNumGEvals(solver, &nge);
	if (flag!=CV_SUCCESS) printf("failed to update SolverStats");
	
	stats.add_counts(nst, nfe, nni, ncfn, netf, njtv, nge);
}

void CVODESolver::print_stats(void) {
    // Kludgy design, move this to somewhere where it will only be called once even if print_stats is called repeatedly
    // Otherwise we'll get double counts of operations
    update_stats(_solver_rupture, _stats_rupture);
    update_stats(_solver_long, _stats_long);
    
	printf("\nSolver Statistics\t\tLong Term\tRupture\n");
	printf("NumSteps\t\t\t\t%-6ld\t\t%-6ld\n", _stats_long._nst, _stats_rupture._nst);
	printf("NumRhsEvals\t\t\t\t%-6ld\t\t%-6ld\n", _stats_long._nfe, _stats_rupture._nfe);
	printf("NumErrTestFails\t\t\t%-6ld\t\t%-6ld\n", _stats_long._netf, _stats_rupture._netf);
	printf("NumNonlinSolvIters\t\t%-6ld\t\t%-6ld\n", _stats_long._nni, _stats_rupture._nni);
	printf("NumNonlinSolvConvFails\t%-6ld\t\t%-6ld\n", _stats_long._ncfn, _stats_rupture._ncfn);
	printf("NumJtimesEvals\t\t\t%-6ld\t\t%-6ld\n", _stats_long._njtv, _stats_rupture._njtv);
	printf("NumRootFindEvals\t\t%-6ld\t\t%-6ld\n", _stats_long._nge, _stats_rupture._nge);
}

int solve_odes(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
	ViCaRS			*sim = static_cast<ViCaRS*>(user_data);
	SimEquations	*eqns = sim->solver()->equations();
	
	// Check if this set of input values is valid, if not return an error
	if (!eqns->values_valid(sim, y)) return 1;
	
	// Solve the equations using the specified method and return the resulting error code
	return eqns->solve_odes(sim, t, y, ydot);
}

int jacobian_times_vector(N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, void *user_data, N_Vector tmp) {
	ViCaRS			*sim = static_cast<ViCaRS*>(user_data);
	SimEquations	*eqns = sim->solver()->equations();
	
	// Solve the equations using the specified method and return the resulting error code
	return eqns->jacobian_times_vector(sim, v, Jv, t, y, fy, tmp);
}

int check_for_mode_change(realtype t, N_Vector y, realtype *gout, void *user_data) {
	ViCaRS			*sim = static_cast<ViCaRS*>(user_data);
	SimEquations	*eqns = sim->solver()->equations();
	
	// Solve the equations using the specified method and return the resulting error code
	return eqns->check_for_mode_change(sim, t, y, gout);
}
