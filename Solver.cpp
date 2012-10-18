#include "Solver.h"
#include "RateState.h"

int RK4Solver::init_solver(ViCaRS *sim) {
	BlockMap::const_iterator	it;
	BlockGID		gid;
    N_Vector        temp_tols;
    
    _num_evals = 0;
    
    k1 = N_VClone(sim->vars());
    k2 = N_VClone(sim->vars());
    k3 = N_VClone(sim->vars());
    k4 = N_VClone(sim->vars());
    in_vals = N_VClone(sim->vars());
    subsum = N_VClone(sim->vars());
    temp_tols = N_VClone(sim->vars());
	
	for (it=sim->begin();it!=sim->end();++it) {
		gid = it->first;
		_eqns->init_block(gid, it->second, sim->vars(), temp_tols);
	}
    
    N_VDestroy_Serial(temp_tols);
    
    return 0;
}

int RK4Solver::advance(ViCaRS *sim, N_Vector vars, realtype target_time, realtype &finish_time, bool &next_solver) {
    realtype        t, dt;
    int             res;
    
    next_solver = 0;
    
    if (!_eqns->values_valid(sim, vars)) return 1;
    
    t = target_time;
    dt = _timestep;
    
    // Calculate k1 = dt * f(y)
    _num_evals++;
    res = _eqns->solve_odes(sim, t, vars, k1);
    if (res) return res;
    N_VScale(dt, k1, k1);
    
    // Calculate k2 = f(y+0.5*k1)
    N_VLinearSum(1, vars, 0.5, k1, in_vals);
    _num_evals++;
    res = _eqns->solve_odes(sim, t+0.5*dt, in_vals, k2);
    if (res) return res;
    N_VScale(dt, k2, k2);
    
    // Calculate k3 = f(y+0.5*k2)
    N_VLinearSum(1, vars, 0.5, k2, in_vals);
    _num_evals++;
    res = _eqns->solve_odes(sim, t+0.5*dt, in_vals, k3);
    if (res) return res;
    N_VScale(dt, k3, k3);
    
    // Calculate k4 = f(y+k3)
    N_VLinearSum(1, vars, 1, k3, in_vals);
    _num_evals++;
    res = _eqns->solve_odes(sim, t+dt, in_vals, k4);
    if (res) return res;
    N_VScale(dt, k4, k4);
    
    // y_{n+1} = y_n + 1/6 * (k1+2k2+2k3+k4)
    //std::cerr << "v: " << NV_Ith_S(vars,1) << " h: " << NV_Ith_S(vars,2) << std::endl;
    //std::cerr << "dv: " << NV_Ith_S(k1,1) << " dh: " << NV_Ith_S(k1,2) << std::endl;
    //std::cerr << "dv: " << NV_Ith_S(k2,1) << " dh: " << NV_Ith_S(k2,2) << std::endl;
    //std::cerr << "dv: " << NV_Ith_S(k3,1) << " dh: " << NV_Ith_S(k3,2) << std::endl;
    //std::cerr << "dv: " << NV_Ith_S(k4,1) << " dh: " << NV_Ith_S(k4,2) << std::endl;
    N_VLinearSum(1, k1, 2, k2, subsum);
    N_VLinearSum(1, subsum, 2, k3, subsum);
    N_VLinearSum(1, subsum, 1, k4, subsum);
    N_VLinearSum(1, vars, 1.0/6.0, subsum, vars);
    
    //std::cerr << "x: " << NV_Ith_S(vars,0) << " v: " << NV_Ith_S(vars,1) << " h: " << NV_Ith_S(vars,2) << std::endl;
    
    finish_time = target_time + dt;
    
    return 0;
}

void RK4Solver::print_stats(void) {
	printf("Number of evaluations: %d\n", _num_evals);
}

CVODESolver::~CVODESolver(void) {
	N_VDestroy_Serial(_abs_tol);
    
	// Free integrator memory
	if (_solver) CVodeFree(&_solver);
}

int CVODESolver::init_solver(ViCaRS *sim) {
	int				flag;
	BlockMap::const_iterator	it;
	BlockGID		gid;
	
	// Initialize variables and tolerances
	_rel_tol = RCONST(1.0e-6);
	_abs_tol = N_VNew_Serial(sim->num_global_blocks()*_eqns->num_equations());
	if (_abs_tol == NULL) return 1;
	
	for (it=sim->begin();it!=sim->end();++it) {
		gid = it->first;
		_eqns->init_block(gid, it->second, sim->vars(), _abs_tol);
	}
    
	// Initialize the solver with a root finder
	flag = init_cvode_solver(&_solver, _rootdir, sim);
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
	flag = CVodeSetMaxNumSteps(solver, 1e5);
	if (flag != CV_SUCCESS) return 1;
    
	// Set the Jacobian x vector function
	if (_eqns->has_jacobian()) {
		flag = CVSpilsSetJacTimesVecFn(solver, jacobian_times_vector);
		if (flag != CV_SUCCESS) return 1;
	}
	
	*created_solver = solver;
	
	return 0;
}

int CVODESolver::advance(ViCaRS *sim, N_Vector vars, realtype target_time, realtype &finish_time, bool &next_solver) {
	int	flag;
	
    next_solver = 0;
    
	flag = CVode(_solver, target_time+_timestep, vars, &finish_time, CV_NORMAL);
	if (flag == CV_SUCCESS) {
		return 0;
	} else if (flag == CV_ROOT_RETURN) {
		// Update the stats for each solver
		update_stats(_solver, _stats);
		
		_eqns->handle_mode_change(sim, finish_time, vars);
		
        next_solver = 1;
        
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
    update_stats(_solver, _stats);
    
	printf("\nSolver Statistics\n");
	printf("NumSteps\t\t\t\t%-6ld\n", _stats._nst);
	printf("NumRhsEvals\t\t\t\t%-6ld\n", _stats._nfe);
	printf("NumErrTestFails\t\t\t%-6ld\n", _stats._netf);
	printf("NumNonlinSolvIters\t\t%-6ld\n", _stats._nni);
	printf("NumNonlinSolvConvFails\t%-6ld\n", _stats._ncfn);
	printf("NumJtimesEvals\t\t\t%-6ld\n", _stats._njtv);
	printf("NumRootFindEvals\t\t%-6ld\n", _stats._nge);
}

int solve_odes(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
	ViCaRS			*sim = static_cast<ViCaRS*>(user_data);
	SimEquations	*eqns = sim->equations();
	
	// Check if this set of input values is valid, if not return an error
	if (!eqns->values_valid(sim, y)) return 1;
	
	// Solve the equations using the specified method and return the resulting error code
	return eqns->solve_odes(sim, t, y, ydot);
}

int jacobian_times_vector(N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, void *user_data, N_Vector tmp) {
	ViCaRS			*sim = static_cast<ViCaRS*>(user_data);
	SimEquations	*eqns = sim->equations();
	
	// Solve the equations using the specified method and return the resulting error code
	return eqns->jacobian_times_vector(sim, v, Jv, t, y, fy, tmp);
}

int check_for_mode_change(realtype t, N_Vector y, realtype *gout, void *user_data) {
	ViCaRS			*sim = static_cast<ViCaRS*>(user_data);
	SimEquations	*eqns = sim->equations();
	
	// Solve the equations using the specified method and return the resulting error code
	return eqns->check_for_mode_change(sim, t, y, gout);
}
