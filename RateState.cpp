#include "RateState.h"

ViCaRS::ViCaRS(unsigned int total_num_blocks) : _num_global_blocks(total_num_blocks), greens_matrix(total_num_blocks, total_num_blocks)
{
    _solver_long = _solver_rupture = _current_solver = NULL;
    _use_slowness_law = true;
	_G = 3.0e10;
}

int ViCaRS::add_block(const BlockGID &id, const BlockData &block_data) {
	_bdata[id] = block_data;
	
	return 0;
}

int ViCaRS::init(EqnSolver *eqns) {
	unsigned int	num_local;
	int				flag, rootdir;
	BlockGID		gid;
	BlockMap::const_iterator	it;
	
	_eqns = eqns;
	_num_equations = _eqns->num_equations();
	
	flag = fill_greens_matrix();
	if (flag != 0) return 1;
	
	num_local = (unsigned int)_bdata.size();
	_vars = N_VNew_Serial(_num_global_blocks*_num_equations);
	if (_vars == NULL) return 1;
	_abs_tol = N_VNew_Serial(_num_global_blocks*_num_equations);
	if (_abs_tol == NULL) return 1;
	
	// Initialize variables and tolerances
	_rel_tol = RCONST(1.0e-7);
	
	for (it=_bdata.begin();it!=_bdata.end();++it) {
		gid = it->first;
		_eqns->init_block(gid, it->second, _vars, _abs_tol);
	}
	
	// Init CVodes structures
	flag = _eqns->init(this);
	if (flag) return flag;
	
	_cur_time = 0;
	
	// Initialize the long term solver
	rootdir = 1;
	flag = init_solver(&_solver_long, rootdir);
	if (flag) return flag;
	
	// And the rupture solver
	rootdir = -1;
	flag = init_solver(&_solver_rupture, rootdir);
	if (flag) return flag;
	
	// Set the current solver to be the long term
	_current_solver = _solver_long;
	
	_in_rupture = false;
	
	return 0;
}

int ViCaRS::init_solver(void **created_solver, int rootdir) {
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
	flag = CVodeInit(solver, solve_odes, _cur_time, _vars);
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
	flag = CVodeSetUserData(solver, this);
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

int ViCaRS::advance(void) {
	int	flag;
	realtype tstep;
	
	if (_in_rupture) tstep = _rupture_timestep;
	else tstep = _long_timestep;
	
	flag = CVode(_current_solver, _cur_time+tstep, _vars, &_cur_time, CV_NORMAL);
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
		flag = CVodeReInit(_solver_long, _cur_time, _vars);
		if (flag != CV_SUCCESS) return 1;
		flag = CVodeReInit(_solver_rupture, _cur_time, _vars);
		if (flag != CV_SUCCESS) return 1;
		_in_rupture = !_in_rupture;
		
		_eqns->handle_mode_change(this, _cur_time, _vars);
		
		return 0;
	} else return flag;
}

void ViCaRS::cleanup(void) {
	// Free y and abstol vectors
	N_VDestroy_Serial(_vars);
	N_VDestroy_Serial(_abs_tol);
	
	// Free integrator memory
	CVodeFree(&_solver_long);
	// Free integrator memory
	CVodeFree(&_solver_rupture);
	
	delete _eqns;
}

realtype ViCaRS::interaction(BlockGID i, BlockGID j) { return greens_matrix.val(i,j); };

void ViCaRS::write_header(FILE *fp) {
	BlockMap::const_iterator	it;
	BlockGID		gid;
	unsigned int		vnum;
	
	fprintf(fp, "t ");
	for (it=_bdata.begin();it!=_bdata.end();++it) {
		gid = it->first;
		for (vnum=0;vnum<_eqns->num_outputs();++vnum) {
			fprintf(fp, "%s%d ", _eqns->var_name(vnum).c_str(), gid);
		}
	}
	fprintf(fp, "\n");
}

void ViCaRS::write_summary_header(FILE *fp) {
	fprintf(fp, "Time\t\tRupture\tV_min(GID)\t\tV_max(GID)\t\n");
}

void ViCaRS::write_summary(FILE *fp) {
	realtype		v, min_v, max_v;
	BlockGID		min_v_gid, max_v_gid;
	BlockMap::const_iterator	it;
	
	min_v = DBL_MAX;
	max_v = -DBL_MAX;
	
	for (it=begin();it!=end();++it) {
		v = 0;//V(it->first);
		if (min_v >= v) {
			min_v = v;
			min_v_gid = it->first;
		}
		if (max_v <= v) {
			max_v = v;
			max_v_gid = it->first;
		}
	}
	fprintf(fp, "%0.2e\t%d\t\t%0.2e(%d)\t%0.2e(%d)\n", _cur_time, _in_rupture, min_v, min_v_gid, max_v, max_v_gid);
}

void ViCaRS::write_cur_data(FILE *fp) {
	BlockMap::const_iterator	it;
	unsigned int				vnum;
	
	fprintf(fp, "%0.7e ", _cur_time);
	for (it=_bdata.begin();it!=_bdata.end();++it) {
		for (vnum=0;vnum<_eqns->num_outputs();++vnum) {
			fprintf(fp, "%14.6e ", _eqns->var_value(this, vnum, it->first, _vars));
		}
	}
	fprintf(fp, "\n");
}

void ViCaRS::update_stats(void *solver, SolverStats &stats) {
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

void ViCaRS::print_stats(void) {
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

int ViCaRS::fill_greens_matrix(void) {
	int i, j, r;
	BlockMap::const_iterator	it, jt;
	
	double G = 3e10;
	double L = 0.1;
	double scale = 1e-9;
	
	for (it=_bdata.begin(); it != _bdata.end(); ++it) {
		i = it->first;
		for (jt=_bdata.begin(); jt != _bdata.end(); ++jt) {
			j = jt->first;
			if (i==j) { greens_matrix.setVal(i,j,scale*-0.53*G/L); }
			else {
				r = abs(i-j);
				greens_matrix.setVal(i,j,scale*0.98*G/(L*pow(((r/L)-1),3)));
			}
		}
	}
	
	return 0;
}

int solve_odes(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
	ViCaRS			*sim = static_cast<ViCaRS*>(user_data);
	EqnSolver		*eqns = sim->get_eqns();
	
	// Check if this set of input values is valid, if not return an error
	if (!eqns->values_valid(sim, y)) return 1;
	
	// Solve the equations using the specified method and return the resulting error code
	return eqns->solve_odes(sim, t, y, ydot);
}

int jacobian_times_vector(N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, void *user_data, N_Vector tmp) {
	ViCaRS			*sim = static_cast<ViCaRS*>(user_data);
	EqnSolver		*eqns = sim->get_eqns();
	
	// Solve the equations using the specified method and return the resulting error code
	return eqns->jacobian_times_vector(sim, v, Jv, t, y, fy, tmp);
}

int check_for_mode_change(realtype t, N_Vector y, realtype *gout, void *user_data) {
	ViCaRS			*sim = static_cast<ViCaRS*>(user_data);
	EqnSolver		*eqns = sim->get_eqns();
	
	// Solve the equations using the specified method and return the resulting error code
	return eqns->check_for_mode_change(sim, t, y, gout);
}

