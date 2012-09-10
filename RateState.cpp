#include "RateState.h"

ViCaRS::ViCaRS(unsigned int total_num_blocks) : _num_global_blocks(total_num_blocks), greens_matrix(total_num_blocks, total_num_blocks)
{
    _solver_long = _solver_rupture = _current_solver = NULL;
    _use_slowness_law = true;
    _use_simple_equations = true;
	_G = 3.0e10;
    //v_ss = 0.01;
    //h_ss = 1/v_ss;
    //v_eq = 0.1;
    //v_min = 1e-9; // values in cvodes vectors cannot be zero, so use this as our 'zero'?
}

int ViCaRS::add_block(const BlockGID &id, const BlockData &block_data) {
	_bdata[id] = block_data;
	
	return 0;
}

int ViCaRS::init(void) {
	unsigned int	num_local;
	int				flag, rootdir;
	BlockGID		gid;
	BlockMap::const_iterator	it;
	
	if (_use_simple_equations) _eqns = new SimpleEqns;
	else _eqns = new OrigEqns;
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
	
	fprintf(fp, "t S ");
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
	
	fprintf(fp, "%0.7e %d ", _cur_time, _in_rupture);
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

/*
 * Functions called by the solver for advanced RS
 */

bool OrigEqns::values_valid(ViCaRS *sim, N_Vector y) {
	int			local_fail;
	BlockMap::const_iterator	it;
	
	if (_use_log_spline) return true;
	
	// Check if any velocity or theta values are below 0
	local_fail = 0;
	for (it=sim->begin();it!=sim->end();++it) {
		if (Vth(y,it->first) <= 0 || Hth(y,it->first) <= 0) {
			local_fail = 1;
			break;
		}
	}
	
	return (!local_fail);
}

std::string OrigEqns::var_name(unsigned int var_num) const {
	switch (var_num) {
		case 0: return "X";
		case 1: return "V";
		case 2: return "H";
		case 3: return "F";
		default: throw std::exception();
	}
}

std::string SimpleEqns::var_name(unsigned int var_num) const {
	switch (var_num) {
		case 0: return "X";
		case 1: return "H";
		case 2: return "V";
		default: throw std::exception();
	}
}

realtype OrigEqns::var_value(ViCaRS *sim, unsigned int var_num, BlockGID gid, N_Vector y) {
	switch (var_num) {
		case 0: return Xth(y, gid);
		case 1: return Vth(y, gid);
		case 2: return Hth(y, gid);
		case 3: return F(sim->param_a(gid), sim->param_a(gid), Vth(y, gid), Hth(y, gid));
		default: throw std::exception();
	}
}

realtype SimpleEqns::var_value(ViCaRS *sim, unsigned int var_num, BlockGID gid, N_Vector y) {
	switch (var_num) {
		case 0:
		case 1:
			return NV_DATA_S(y)[gid*num_equations()+var_num];
		case 2:
			return 0;	// Calculate V
		default: throw std::exception();
	}
}

void OrigEqns::init_block(BlockGID gid, const BlockData &block, N_Vector vars, N_Vector tols) {
	Xth(vars, gid) = block._init_x;
	Vth(vars, gid) = block._init_v;
	Hth(vars, gid) = block._init_h;
	
	Xth(tols, gid) = block._tol_x;
	Vth(tols, gid) = block._tol_v;
	Hth(tols, gid) = block._tol_h;
}

void SimpleEqns::init_block(BlockGID gid, const BlockData &block, N_Vector vars, N_Vector tols) {
	Xth(vars, gid) = block._init_x;
	Hth(vars, gid) = block._init_h;
	
	Xth(tols, gid) = block._tol_x;
	Hth(tols, gid) = block._tol_h;
}

int OrigEqns::init(ViCaRS *sim) {
	return 0;
}

int SimpleEqns::init(ViCaRS *sim) {
	BlockMap::const_iterator	it,jt;
	realtype					v_ss, sum_load, W, beta, delta_tau;
	
	// Start all blocks in phase 0
	for (it=sim->begin();it!=sim->end();++it) {
		phase[it->first] = 0;
	}
	
	mu_0 = 0.5;
	A = 0.005;
	B = 0.015;
	D_c = 0.01*1e-3;
	sigma_i = 1;
	_v_star = 1.0/(1.0e5*365.25*86400);	// 1 meter/1e5 years in m/s
	_theta_star = D_c/_v_star;
	
	v_ss = 1.0/(1.0e2*365.25*86400);	// 1 cm/year in m/s
	W = 1000;
	
	beta = 5000;		// meters/second
	delta_tau = 1000;
	
	for (it=sim->begin();it!=sim->end();++it) {
		_ss_stress[it->first] = sigma_i*(mu_0+(A-B)*log(v_ss/_v_star));
		_v_eq[it->first] = 2*beta*delta_tau/sim->G();
		_base_stress[it->first] = _elem_stress[it->first] = 0;
	}
	
	for (it=sim->begin();it!=sim->end();++it) {
		sum_load = 0;
		for (jt=sim->begin();jt!=sim->end();++jt) {
			sum_load += v_ss*sim->interaction(it->first, jt->first);
		}
		sum_load += v_ss*sim->G()/W;
		_stress_loading[it->first] = sum_load;
	}
	
	return 0;
}

int OrigEqns::solve_odes(ViCaRS *sim, realtype t, N_Vector y, N_Vector ydot) {
	realtype		friction_force, spring_force, x, v, h, interact;
	unsigned int	i_gid, j_gid;
	BlockMap::const_iterator	it,jt;
	
	// Share X values over all processors to enable correct interaction between blocks
	// TODO: only share X values with processors that need them
	
	for (it=sim->begin();it!=sim->end();++it) {
		i_gid = it->first;
		x = Xth(y,i_gid);
		v = Vth(y,i_gid);
		h = Hth(y,i_gid);
		
		friction_force = sim->param_k(i_gid)*F(sim->param_a(i_gid),sim->param_b(i_gid),v,h);
		
		// calculate interaction force
		spring_force = 0;
		for (jt=sim->begin();jt!=sim->end();++jt) {
			j_gid = jt->first;
			interact = sim->interaction(i_gid,j_gid);
			if (interact > 0) {
				spring_force += interact*(x-Xth(y,j_gid));
			}
		}
		spring_force = 0;
		
		Xth(ydot,i_gid) = v;
		Vth(ydot,i_gid) = (t-x-friction_force-spring_force)/sim->param_r(i_gid);
        if (sim->use_slowness_law()) Hth(ydot,i_gid) = RCONST(1) - h*v;
        else Hth(ydot,i_gid) = -h*v*log_func(h*v);
	}
	
	return 0;
}

// Function to determine when the simulation moves into rupture or long term mode
// TODO: parallelize this function
int OrigEqns::check_for_mode_change(ViCaRS *sim, realtype t, N_Vector y, realtype *gout) {
	BlockMap::const_iterator	it;
	realtype			local_max = -DBL_MAX;
	
	for (it=sim->begin();it!=sim->end();++it) local_max = fmax(local_max, Vth(y,it->first));
	
	gout[0] = local_max - sim->rupture_threshold();
	
	return 0;
}

// Computes Jv = Jacobian times v in order to improve convergence for Krylov subspace solver
int OrigEqns::jacobian_times_vector(ViCaRS *sim, N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, N_Vector tmp) {
	BlockMap::const_iterator	it;
	BlockGID	gid;
	
	for (it=sim->begin();it!=sim->end();++it) {
		gid = it->first;
		Xth(Jv,gid) =	1 * Vth(v,gid);
		Vth(Jv,gid) =	(-RCONST(1)/sim->param_r(gid)) * Xth(v,gid)
		+ -(sim->param_k(gid)*sim->param_a(gid)*log_deriv((Vth(y,gid)))/(sim->param_r(gid))) * Vth(v,gid)
		+ -(sim->param_k(gid)*sim->param_b(gid)/(sim->param_r(gid)*Hth(y,gid))) * Hth(v,gid);
        
        if (sim->use_slowness_law()) {
            Hth(Jv,gid) =	-Hth(y,gid) * Vth(v,gid)
			+ -Vth(y,gid) * Hth(v,gid);
        } else {
            Hth(Jv,gid) = -Hth(y,gid)*(log_func(Hth(y,gid)*Vth(y,gid))
									   + Hth(y,gid)*Vth(y,gid) * log_deriv(Hth(y,gid)*Vth(y,gid))) * Vth(v,gid)
			- Vth(y,gid)*(log_func(Hth(y,gid)*Vth(y,gid))
						  + Hth(y,gid)*Vth(y,gid)*log_deriv(Hth(y,gid)*Vth(y,gid))) * Hth(v,gid);
        }
	}
	
	return 0;
}

// Simple eqns Jacobian
// Term		Phase 0		Phase 1		Phase 2
// dX/dX =	0			0			0
// dX/dH =	0
// dH/dX = 
int SimpleEqns::solve_odes(ViCaRS *sim, realtype t, N_Vector y, N_Vector ydot) {
	BlockMap::const_iterator	it, j;
	int             phase_num;
	BlockGID        gid;
	realtype        x, h, v_0=0.001, tau_dot, H, K;
	
	// Calculate velocities and stresses on blocks
	for (it=sim->begin();it!=sim->end();++it) {
		gid = it->first;
		switch (phase[gid]) {
			case 0:
				_vel[gid] = 0;
				_elem_stress[gid] = _base_stress[gid] + (t-_start_time[gid])*_stress_loading[gid];
				break;
			case 1:
				K = -sim->interaction(gid, gid)+1;	// K_T = 1 arbitrarily, need to change this
				H = B/D_c-K/sigma_i;
				tau_dot = _stress_loading[gid];
				_vel[gid] = 1.0/((1.0/v_0+H*sigma_i/tau_dot)*exp(-tau_dot*(t-_start_time[gid])/(A*sigma_i)) - H*sigma_i/tau_dot);
				_elem_stress[gid] = sigma_i*(mu_0+A*log(_vel[gid]/_v_star)+B*log(Hth(y,gid)/_theta_star));
				break;
			case 2:
				_vel[gid] = _v_eq[gid];
				_elem_stress[gid] = sigma_i*(mu_0+A*log(_vel[gid]/_v_star)+B*log(Hth(y,gid)/_theta_star));
				break;
		}
	}
	
	for (it=sim->begin();it!=sim->end();++it) {
		gid = it->first;
		
		phase_num = phase[gid];
		
		x = Xth(y,gid);
		h = Hth(y,gid);
		
		if (phase_num == 0) {
			Hth(ydot,gid) = 1;
		} else if (phase_num == 1) {
			Hth(ydot,gid) = 1.0-h*_vel[gid]/D_c;
		} else if (phase_num == 2) {
			Hth(ydot,gid) = 1.0-h*_vel[gid]/D_c;
		} else return -1;
		
		Xth(ydot,gid) = _vel[gid];
	}
	
	return 0;
}

// Function to determine when the simulation moves into rupture or long term mode
int SimpleEqns::check_for_mode_change(ViCaRS *sim, realtype t, N_Vector y, realtype *gout) {
	BlockMap::const_iterator it;
	realtype	max_val;
	BlockGID gid;
	
	max_val = -DBL_MAX;
	for (it=sim->begin();it!=sim->end();++it) {
		gid = it->first;
		
		if (phase[gid] == 0) max_val = fmax(max_val, _elem_stress[gid] - _ss_stress[gid]);
		else if (phase[gid] == 1) max_val = fmax(max_val, _vel[gid] - _v_eq[gid]);
		else if (phase[gid] == 2) max_val = fmax(max_val, _ss_stress[gid] - _elem_stress[gid]);
		else return 1;
	}
	
	gout[0] = max_val;
	
	return 0;
}

void OrigEqns::handle_mode_change(ViCaRS *sim, realtype t, N_Vector y) {
}

void SimpleEqns::handle_mode_change(ViCaRS *sim, realtype t, N_Vector y) {
	BlockMap::const_iterator it;
	BlockGID gid;
	
	for (it=sim->begin();it!=sim->end();++it) {
		gid = it->first;
		
		if (phase[gid] == 0 && _elem_stress[gid] >= _ss_stress[gid]) {
			phase[gid] = 1;
		} else if (phase[gid] == 1 && _vel[gid] >= _v_eq[gid]) {
			phase[gid] = 2;
		} else if (phase[gid] == 2 && _ss_stress[gid] >= _elem_stress[gid]) {
			phase[gid] = 0;
		}
	}
}

