#include "RateState.h"

ViCaRS::ViCaRS(unsigned int total_num_blocks) :_num_global_blocks(total_num_blocks), _num_equations(NEQ), greens_matrix(total_num_blocks, total_num_blocks), log_approx(-40, 2, 1e12, 15, 100, 1e5) {
    _solver_long = _solver_rupture = _current_solver = NULL;
    _use_slowness_law = true;
    _use_log_spline = true;
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
	int				flag, rootdir[1];
	realtype		*toldata;
	GlobalLocalMap::const_iterator	it;

  flag = fill_greens_matrix();
  if (flag != 0) return 1;
	
	num_local = (unsigned int)_global_local_map.size();
	_vars = N_VNew_Parallel(MPI_COMM_WORLD, num_local*_num_equations, _num_global_blocks*_num_equations);
	if (_vars == NULL) return 1;
	_abs_tol = N_VNew_Parallel(MPI_COMM_WORLD, num_local*_num_equations, _num_global_blocks*_num_equations);
	if (_abs_tol == NULL) return 1;
	
	// Initialize variables and tolerances
	_rel_tol = RCONST(1.0e-10);
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
	
	/* Call CVodeCreate to create the solver memory and specify the 
	 * Backward Differentiation Formula and the use of a Newton iteration */
	_solver_long = CVodeCreate(CV_BDF, CV_NEWTON);
	if (_solver_long == NULL) return 1;
	_solver_rupture = CVodeCreate(CV_BDF, CV_NEWTON);
	if (_solver_rupture == NULL) return 1;
    
	// Turn off error messages
	//flag = CVodeSetErrFile(_solver_long, NULL);
	//if (flag != CV_SUCCESS) return 1;
	//flag = CVodeSetErrFile(_solver_rupture, NULL);
	//if (flag != CV_SUCCESS) return 1;
	
	// Call CVodeInit to initialize the integrator memory and specify the
	// user's right hand side function in y'=f(t,y), the inital time, and
	// the initial dependent variable vector.
	_cur_time = 0;
	flag = CVodeInit(_solver_long, func, _cur_time, _vars);
	if (flag != CV_SUCCESS) return 1;
	flag = CVodeInit(_solver_rupture, func, _cur_time, _vars);
	if (flag != CV_SUCCESS) return 1;
	
	// Call CVodeSVtolerances to specify the scalar relative tolerance
	// and vector absolute tolerances
	flag = CVodeSVtolerances(_solver_long, _rel_tol, _abs_tol);
	if (flag != CV_SUCCESS) return 1;
	flag = CVodeSVtolerances(_solver_rupture, _rel_tol, _abs_tol);
	if (flag != CV_SUCCESS) return 1;
	
	// Set the root finding function, going above the rupture limit for long term solver
	// and going below the limit for rupture solver.
	flag = CVodeRootInit(_solver_long, 1, rupture_test);
	if (flag != CV_SUCCESS) return 1;
	rootdir[0] = 1;
	flag = CVodeSetRootDirection(_solver_long, rootdir);
	if (flag != CV_SUCCESS) return 1;
	
	flag = CVodeRootInit(_solver_rupture, 1, rupture_test);
	if (flag != CV_SUCCESS) return 1;
	rootdir[0] = -1;
	flag = CVodeSetRootDirection(_solver_rupture, rootdir);
	if (flag != CV_SUCCESS) return 1;
	
	// Call CVSpbcg to specify the CVSPBCG scaled preconditioned Bi-CGSTab iterative solver
	// Currently not using any preconditioner
	flag = CVSpbcg(_solver_long, PREC_NONE, 0);
	if (flag != CV_SUCCESS) return 1;
	flag = CVSpbcg(_solver_rupture, PREC_NONE, 0);
	if (flag != CV_SUCCESS) return 1;
	
	// Set the user data
	flag = CVodeSetUserData(_solver_long, this);
	if (flag != CV_SUCCESS) return 1;
	flag = CVodeSetUserData(_solver_rupture, this);
	if (flag != CV_SUCCESS) return 1;
	
	// Set maximum number of steps per solver iteration
	// This should be higher if the time step is less frequent or tolerance is lowered
	flag = CVodeSetMaxNumSteps(_solver_long, 1e7);
	if (flag != CV_SUCCESS) return 1;
	flag = CVodeSetMaxNumSteps(_solver_rupture,1e7);
	if (flag != CV_SUCCESS) return 1;

	// Set the Jacobian x vector function
	flag = CVSpilsSetJacTimesVecFn(_solver_long, jacobian_times_vector);
	if (flag != CV_SUCCESS) return 1;
	flag = CVSpilsSetJacTimesVecFn(_solver_rupture, jacobian_times_vector);
	if (flag != CV_SUCCESS) return 1;
	
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
    // Set the current solver to be the long term
    _current_solver = _solver_long;
    
    _in_rupture = false;

	return 0;
}

int ViCaRS::advance(void) {
	int				flag;
	realtype		tstep;
    
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
	N_VDestroy_Parallel(_abs_tol);
	
	// Free integrator memory
	CVodeFree(&_solver_long);
	// Free integrator memory
	CVodeFree(&_solver_rupture);
}

realtype ViCaRS::interaction(BlockGID a, BlockGID b) {
  return greens_matrix.val(a,b);
}

void ViCaRS::write_header(FILE *fp) {
	GlobalLocalMap::const_iterator	it;
	BlockGID		gid;
	fprintf(fp, "t S ");
	for (it=_global_local_map.begin();it!=_global_local_map.end();++it) {
		gid = it->first;
		fprintf(fp, "x%d v%d h%d F%d ", gid, gid, gid, gid);
	}
	fprintf(fp, "\n");
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
	flag = CVodeGetNumRhsEvals(solver, &nfe);
	flag = CVodeGetNumErrTestFails(solver, &netf);
	flag = CVodeGetNumNonlinSolvIters(solver, &nni);
	flag = CVodeGetNumNonlinSolvConvFails(solver, &ncfn);
	flag = CVSpilsGetNumJtimesEvals(solver, &njtv);
	flag = CVodeGetNumGEvals(solver, &nge);
    
    stats.add_counts(nst, nfe, nni, ncfn, netf, njtv, nge);
}

void ViCaRS::print_stats(void) {
    // Kludgy design, move this to somewhere where it will only be called once even if print_stats is called repeatedly
    // Otherwise we'll get double counts of operations
    update_stats(_solver_rupture, _stats_rupture);
    update_stats(_solver_long, _stats_long);
    
	printf("\nSolver Statistics\t\tLong Term\tRupture\n");
	printf("NumSteps\t\t\t%-6ld\t\t%-6ld\n", _stats_long._nst, _stats_rupture._nst);
	printf("NumRhsEvals\t\t\t%-6ld\t\t%-6ld\n", _stats_long._nfe, _stats_rupture._nfe);
	printf("NumErrTestFails\t\t\t%-6ld\t\t%-6ld\n", _stats_long._netf, _stats_rupture._netf);
	printf("NumNonlinSolvIters\t\t%-6ld\t\t%-6ld\n", _stats_long._nni, _stats_rupture._nni);
	printf("NumNonlinSolvConvFails\t\t%-6ld\t\t%-6ld\n", _stats_long._ncfn, _stats_rupture._ncfn);
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
    i = it->second;
    for (jt=_global_local_map.begin(); jt != _global_local_map.end(); ++jt) {
      j = jt->second;
      if (i==j) { greens_matrix.setVal(i,j,scale*-0.53*G/L); }
      else {
        r = abs(i-j);
        greens_matrix.setVal(i,j,scale*0.98*G/(L*pow(((r/L)-1),3)));
      }
    }
  }

  return 0;
}

/*
 * Functions called by the solver for advanced RS
 */

int func(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
	realtype		friction_force, spring_force, x, v, h, interact;
	int				local_fail, global_fail;
	unsigned int	lid, gid;
	BlockGID		other_block;
	ViCaRS			*sim = (ViCaRS*)(user_data);
	GlobalLocalMap::const_iterator	it;
	realtype		*global_x;
	
	// Check if any velocity or theta values are below 0
	local_fail = 0;
    if (sim->use_log_spline()) {
        for (it=sim->begin();it!=sim->end();++it) {
            if (Vth(y,it->second) <= 0 || Hth(y,it->second) <= 0) {
                local_fail = 1;
                break;
            }
        }
    }

	// Communicate with other processes to indicate whether or not to continue
	// TODO: error check
	MPI_Allreduce(&local_fail, &global_fail, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
	if (global_fail) return 1;
	
	// Share X values over all processors to enable correct interaction between blocks
	// TODO: only share X values with processors that need them
	global_x = new realtype[sim->num_global_blocks()];
	
	for (it=sim->begin();it!=sim->end();++it) {
		gid = it->first;
		lid = it->second;
		x = Xth(y,lid);
		v = Vth(y,lid);
		h = Hth(y,lid);
		
		friction_force = sim->param_k(gid)*sim->F(gid,v,h);
		spring_force = 0;
		for (other_block=0;other_block<sim->num_global_blocks();++other_block) {
			interact = sim->interaction(gid, other_block);
			if (interact > 0) {
				spring_force += interact*(x-Xth(y,sim->global_to_local(other_block)));
			}
		}
		
		Xth(ydot,lid) = v;
		Vth(ydot,lid) = (t-x-friction_force-spring_force)/sim->param_r(gid);
        if (sim->use_slowness_law()) Hth(ydot,lid) = RCONST(1) - h*v;
        else Hth(ydot,lid) = -h*v*sim->log_func(h*v);
	}
	
	delete global_x;
	
	return 0;
}

// Function to determine when the simulation moves into rupture or long term mode
// TODO: parallelize this function
int rupture_test(realtype t, N_Vector y, realtype *gout, void *user_data) {
	GlobalLocalMap::const_iterator	it;
	ViCaRS			*sim = (ViCaRS*)(user_data);
	double			local_max = -DBL_MAX;
	
	for (it=sim->begin();it!=sim->end();++it) local_max = fmax(local_max, Vth(y,it->second));
	
	gout[0] = local_max - sim->rupture_threshold();
	
	return 0;
}

// Computes Jv = Jacobian times v in order to improve convergence for Krylov subspace solver
int jacobian_times_vector(N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, void *user_data, N_Vector tmp) {
	ViCaRS			*sim = (ViCaRS*)(user_data);
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

/*
 * Functions called by the solver for simple RS
 */

int simple_equations(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
    return 0;
}

// Function to determine when the simulation moves into rupture or long term mode
// TODO: parallelize this function
int simple_rupture_test(realtype t, N_Vector y, realtype *gout, void *user_data) {
    return 0;
}

// Computes Jv = Jacobian times v in order to improve convergence for Krylov subspace solver
int simple_jacobian_times_vector(N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, void *user_data, N_Vector tmp) {
    return 0;
}
