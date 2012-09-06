#include "RateState.h"

ViCaRS::ViCaRS(unsigned int total_num_blocks) : _num_global_blocks(total_num_blocks), _num_equations(NEQ), log_approx(-40, 2, 1e12, 1e5, 15), greens_matrix(total_num_blocks, total_num_blocks)
{
    _solver_long = _solver_rupture = _current_solver = NULL;
    _use_slowness_law = true;
    _use_log_spline = false;
    _use_simple_equations = false;
    v_ss = 0.01;
    h_ss = 1/v_ss;
    v_eq = 0.1;
    v_min = 1e-9; // values in cvodes vectors cannot be zero, so use this as our 'zero'?
}

int ViCaRS::add_local_block(const BlockData &block_data) {
	BlockGID		gid;
	
	gid = block_data._gid;
	
	if (gid >= _num_global_blocks) return -1;
	
	if (_global_local_map.count(gid) > 0) {
		_bdata.at(_global_local_map[gid]) = block_data;
	} else {
		_global_local_map[gid] = (unsigned int)_bdata.size();
		_bdata.push_back(block_data);
	}
	
	return 0;
}

int ViCaRS::init(void) {
	unsigned int	num_local, lid;
	int				flag, rootdir;
	realtype		*toldata;
	GlobalLocalMap::const_iterator	it;
	
	flag = fill_greens_matrix();
	if (flag != 0) return 1;
	
	num_local = (unsigned int)_global_local_map.size();
	_vars = N_VNew_Parallel(MPI_COMM_WORLD, num_local*_num_equations, _num_global_blocks*_num_equations);
	if (_vars == NULL) return 1;
	_stress = N_VNew_Parallel(MPI_COMM_WORLD, num_local, _num_global_blocks);
	if (_stress == NULL) return 1;
	_abs_tol = N_VNew_Parallel(MPI_COMM_WORLD, num_local*_num_equations, _num_global_blocks*_num_equations);
	if (_abs_tol == NULL) return 1;
	
	// Initialize variables and tolerances
	_rel_tol = RCONST(1.0e-7);
	toldata = NV_DATA_P(_abs_tol);
	
	for (it=_global_local_map.begin();it!=_global_local_map.end();++it) {
		lid = it->second;
		X(it->first) = _bdata[lid]._init_x;
		V(it->first) = _bdata[lid]._init_v;
		H(it->first) = _bdata[lid]._init_h;
		
		toldata[_num_equations*lid+EQ_X] = _bdata[lid]._tol_x;
		toldata[_num_equations*lid+EQ_V] = _bdata[lid]._tol_v;
		toldata[_num_equations*lid+EQ_H] = _bdata[lid]._tol_h;
	}
	
	// Init CVodes structures
	if (_use_simple_equations) _eqns = new SimpleEqns;
	else _eqns = new OrigEqns;
	
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
	flag = CVodeRootInit(solver, 1, check_for_rupture);
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

/*
int ViCaRS::advance_simple(void) {
	int flag, block_phase;
	realtype tstep;
	unsigned int lid;
	BlockGID gid;
	
	// phase changed occured during last time-step
	else if (flag == CV_ROOT_RETURN) {
		update_stats(_solver_long, _stats_long);
		update_stats(_solver_rupture, _stats_rupture);
		
		CVodeGetRootInfo(_current_solver, _roots);
		
		_in_rupture = false;
		
		GlobalLocalMap::const_iterator	it;
		for (it=_global_local_map.begin();it!=_global_local_map.end();++it) {
			gid = it->first;
			lid = it->second;
			
			// if block phase change detected by cvodes root-finder
			if (_roots[gid] != 0) {
				
				block_phase = phase[gid];
				
				// going into rupture mode: 0 -> 1 or 1 -> 2
				if (block_phase == 0 || block_phase == 1) _in_rupture = true;
				
				// change block phase
				if (block_phase == 0) {
					phase[gid] = 1;
					Vth(_vars,lid) = v_ss;
					Hth(_vars,lid) = h_ss;
				}
				else if (block_phase == 1) {
					phase[gid] = 2;
					Vth(_vars,lid) = v_eq;
				}
				else if (block_phase == 2) {
					phase[gid] = 0;
					Vth(_vars,lid) = v_min;
					Hth(_vars,lid) = 1;
				}
				else return -1;
			}
		}
		
		return 0;
	}
	
	else return flag;
}*/

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
		
		return 0;
	} else return flag;
}

void ViCaRS::cleanup(void) {
	// Free y and abstol vectors
	N_VDestroy_Parallel(_vars);
	N_VDestroy_Parallel(_stress);
	N_VDestroy_Parallel(_abs_tol);
	
	// Free integrator memory
	CVodeFree(&_solver_long);
	// Free integrator memory
	CVodeFree(&_solver_rupture);
	
	delete _eqns;
}

realtype ViCaRS::interaction(BlockGID i, BlockGID j) { return greens_matrix.val(i,j); };

void ViCaRS::write_header(FILE *fp) {
	GlobalLocalMap::const_iterator	it;
	BlockGID		gid;
	int world_size, rank;
	
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	if (rank == 0) {
		fprintf(fp, "t S ");
		for (it=_global_local_map.begin();it!=_global_local_map.end();++it) {
			gid = it->first;
			fprintf(fp, "x%d v%d h%d F%d ", gid, gid, gid, gid);
		}
		fprintf(fp, "\n");
	}
}

void ViCaRS::write_summary_header(FILE *fp) {
	fprintf(fp, "Time\t\tRupture\tV_min(GID)\t\tV_max(GID)\t\n");
}

void ViCaRS::write_summary(FILE *fp) {
	realtype		v, min_v, max_v;
	BlockGID		min_v_gid, max_v_gid;
	GlobalLocalMap::const_iterator	it;
	
	min_v = DBL_MAX;
	max_v = -DBL_MAX;
	
	for (it=begin();it!=end();++it) {
		v = V(it->first);
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
	realtype			xi, vi, hi;
	GlobalLocalMap::const_iterator	it;
	
	fprintf(fp, "%0.7e %d ", _cur_time, _in_rupture);
	for (it=_global_local_map.begin();it!=_global_local_map.end();++it) {
		xi = X(it->first);
		vi = V(it->first);
		hi = H(it->first);
		fprintf(fp, "%14.6e %14.6e %14.6e %14.6e ", xi, vi, hi, F(it->first, vi, hi));
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
	GlobalLocalMap::const_iterator	it, jt;
	
	double G = 3e10;
	double L = 0.1;
	double scale = 1e-9;
	
	for (it=_global_local_map.begin(); it != _global_local_map.end(); ++it) {
		i = it->first;
		for (jt=_global_local_map.begin(); jt != _global_local_map.end(); ++jt) {
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

int check_for_rupture(realtype t, N_Vector y, realtype *gout, void *user_data) {
	ViCaRS			*sim = static_cast<ViCaRS*>(user_data);
	EqnSolver		*eqns = sim->get_eqns();
	
	// Solve the equations using the specified method and return the resulting error code
	return eqns->check_for_rupture(sim, t, y, gout);
}

/*
 * Functions called by the solver for advanced RS
 */

bool EqnSolver::values_valid(ViCaRS *sim, N_Vector y) {
	int			local_fail, global_fail, err;
	GlobalLocalMap::const_iterator	it;
	
	if (sim->use_log_spline()) return true;
	
	// Check if any velocity or theta values are below 0
	local_fail = 0;
	for (it=sim->begin();it!=sim->end();++it) {
		if (Vth(y,it->second) <= 0 || Hth(y,it->second) <= 0) {
			local_fail = 1;
			break;
		}
	}
	
	// Communicate with other processes to indicate whether or not to continue
	// TODO: error check
	err = MPI_Allreduce(&local_fail, &global_fail, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
	return (!global_fail);
}

int OrigEqns::solve_odes(ViCaRS *sim, realtype t, N_Vector y, N_Vector ydot) {
	realtype		friction_force, spring_force, x, v, h, interact;
	unsigned int	i_lid, j_lid, i_gid;
	GlobalLocalMap::const_iterator	it,jt;
	
	// Share X values over all processors to enable correct interaction between blocks
	// TODO: only share X values with processors that need them
	
	for (it=sim->begin();it!=sim->end();++it) {
		i_gid = it->first;
		i_lid = it->second;
		x = Xth(y,i_lid);
		v = Vth(y,i_lid);
		h = Hth(y,i_lid);
		
		friction_force = sim->param_k(i_gid)*sim->F(i_gid,v,h);
		
		// calculate interaction force
		spring_force = 0;
		for (jt=sim->begin();jt!=sim->end();++jt) {
			j_lid = jt->second;
			interact = sim->interaction(i_lid,j_lid);
			if (interact > 0) {
				spring_force += interact*(x-Xth(y,j_lid));
			}
		}
		spring_force = 0;
		
		Xth(ydot,i_lid) = v;
		Vth(ydot,i_lid) = (t-x-friction_force-spring_force)/sim->param_r(i_gid);
        if (sim->use_slowness_law()) Hth(ydot,i_lid) = RCONST(1) - h*v;
        else Hth(ydot,i_lid) = -h*v*sim->log_func(h*v);
	}
	
	return 0;
}

// Function to determine when the simulation moves into rupture or long term mode
// TODO: parallelize this function
int OrigEqns::check_for_rupture(ViCaRS *sim, realtype t, N_Vector y, realtype *gout) {
	GlobalLocalMap::const_iterator	it;
	double			local_max = -DBL_MAX;
	realtype		local_gout;
	
	for (it=sim->begin();it!=sim->end();++it) local_max = fmax(local_max, Vth(y,it->second));
	
	local_gout = local_max - sim->rupture_threshold();
	
	MPI_Allreduce(&local_gout, gout, 1, PVEC_REAL_MPI_TYPE, MPI_MIN, MPI_COMM_WORLD);
	
	return 0;
}

// Computes Jv = Jacobian times v in order to improve convergence for Krylov subspace solver
int OrigEqns::jacobian_times_vector(ViCaRS *sim, N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, N_Vector tmp) {
	GlobalLocalMap::const_iterator	it;
	unsigned int	lid, gid;
	
	for (it=sim->begin();it!=sim->end();++it) {
		gid = it->first;
		lid = it->second;
		Xth(Jv,lid) =	1 * Vth(v,lid);
		Vth(Jv,lid) =	(-RCONST(1)/sim->param_r(gid)) * Xth(v,lid)
		+ -(sim->param_k(gid)*sim->param_a(gid)*sim->log_deriv((Vth(y,lid)))/(sim->param_r(gid))) * Vth(v,lid)
		+ -(sim->param_k(gid)*sim->param_b(gid)/(sim->param_r(gid)*Hth(y,lid))) * Hth(v,lid);
        
        if (sim->use_slowness_law()) {
            Hth(Jv,lid) =	-Hth(y,lid) * Vth(v,lid)
			+ -Vth(y,lid) * Hth(v,lid);
        } else {
            Hth(Jv,lid) = -Hth(y,lid)*(sim->log_func(Hth(y,lid)*Vth(y,lid))
									   + Hth(y,lid)*Vth(y,lid) * sim->log_deriv(Hth(y,lid)*Vth(y,lid))) * Vth(v,lid)
			- Vth(y,lid)*(sim->log_func(Hth(y,lid)*Vth(y,lid))
						  + Hth(y,lid)*Vth(y,lid)*sim->log_deriv(Hth(y,lid)*Vth(y,lid))) * Hth(v,lid);
        }
		//std::cerr << lid << " " << gid << " " << Xth(Jv,lid) << " " << Vth(Jv,lid) << " " << Hth(Jv,lid) << std::endl;
	}
	
	return 0;
}

int SimpleEqns::solve_odes(ViCaRS *sim, realtype t, N_Vector y, N_Vector ydot) {
	GlobalLocalMap::const_iterator	it, j;
	unsigned int	  lid;
	int             phase_num;
	BlockGID        gid;
	realtype        x, v, h, sigma, mu;
	
	for (it=sim->begin();it!=sim->end();++it) {
		gid = it->first;
		lid = it->second;
		
		phase_num = phase[gid];
		x = Xth(y,lid);
		v = Vth(y,lid);
		h = Hth(y,lid);
		
		/*for (it=sim->begin();it!=sim->end();++it) {
		 mu = 1;
		 sigma = 1;
		 NV_Ith_P(sim->_stress,lid) = sigma*(mu
		 +(phase != 0 ? sim->param_a(gid)*log(v) : 0)
		 +sim->param_b(gid)*log(h));
		 }*/
		
		Xth(ydot,lid) = v;
		
		if (phase_num == 0) {
			Vth(ydot,lid) = 0;
			Hth(ydot,lid) = 1;
		} else if (phase_num == 1) {
			Vth(ydot,lid) = t-x-sim->param_k(gid)*sim->F(gid, v, h);
			Hth(ydot,lid) = 1.0-h*v;
		} else if (phase_num == 2) {
			Vth(ydot,lid) = 0;
			Hth(ydot,lid) = 1.0-h*v;
		} else return -1;
	}
	
	return 0;
}

// Function to determine when the simulation moves into rupture or long term mode
int SimpleEqns::check_for_rupture(ViCaRS *sim, realtype t, N_Vector y, realtype *gout) {
	GlobalLocalMap::const_iterator it;
	int phase_num;
	unsigned int lid;
	BlockGID gid;
	
	for (it=sim->begin();it!=sim->end();++it) {
		gid = it->first;
		lid = it->second;
		
		phase_num = phase[gid];
		if (phase_num == 0) gout[gid] = Hth(y,lid) - sim->h_ss;
		else if (phase_num == 1) gout[gid] = Vth(y,lid) - sim->v_eq;
		else if (phase_num == 2) gout[gid] = sim->h_ss - Hth(y,lid);
		else return 1;
	}
	
	return 0;
}
