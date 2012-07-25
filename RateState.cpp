#include "RateState.h"

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
	int				flag;
	realtype		*toldata;
	GlobalLocalMap::const_iterator	it;
	
	num_local = (unsigned int)_global_local_map.size();
	_vars = N_VNew_Parallel(MPI_COMM_WORLD, num_local*_num_equations, _num_global_blocks*_num_equations);
	if (_vars == NULL) return 1;
	_abs_tol = N_VNew_Parallel(MPI_COMM_WORLD, num_local*_num_equations, _num_global_blocks*_num_equations);
	if (_abs_tol == NULL) return 1;
	
	// Initialize variables and tolerances
	_rel_tol = RCONST(1.0e-8);
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
	_solver = CVodeCreate(CV_BDF, CV_NEWTON);
	if (_solver == NULL) return 1;
	
	// Turn off error messages
	//CVodeSetErrFile(_solver, NULL);
	
	// Call CVodeInit to initialize the integrator memory and specify the
	// user's right hand side function in y'=f(t,y), the inital time, and
	// the initial dependent variable vector.
	_cur_time = 0;
	flag = CVodeInit(_solver, func, _cur_time, _vars);
	if (flag != CV_SUCCESS) return 1;
	
	// Call CVodeSVtolerances to specify the scalar relative tolerance
	// and vector absolute tolerances
	flag = CVodeSVtolerances(_solver, _rel_tol, _abs_tol);
	if (flag != CV_SUCCESS) return 1;
	
	// Set the root finding function
	//flag = CVodeRootInit(cvode, params.num_blocks(), vel_switch_finder);
	//if (check_flag(&flag, "CVodeRootInit", 1)) return(1);
	
	// Call CVSpbcg to specify the CVSPBCG scaled preconditioned Bi-CGSTab iterative solver
	// Currently not using any preconditioner
	flag = CVSpbcg(_solver, PREC_NONE, 0);
	if (flag != CV_SUCCESS) return 1;
	
	// Set the user data
	flag = CVodeSetUserData(_solver, this);
	if (flag != CV_SUCCESS) return 1;
	
	// Set maximum number of steps per solver iteration
	// This should be higher if the time step is less frequent or tolerance is lowered
	flag = CVodeSetMaxNumSteps(_solver, 1000000);
	if (flag != CV_SUCCESS) return 1;

	// Set the Jacobian x vector function
	//flag = CVSpilsSetJacTimesVecFn(_solver, jacobian_times_vector);
	if (flag != CV_SUCCESS) return 1;
	
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	num_send = new int[world_size];
	num_recv = new int[world_size];
	send_offset = new int[world_size];
	recv_offset = new int[world_size];

	return 0;
}

int ViCaRS::advance(realtype next_time) {
	int				flag;
	
	flag = CVode(_solver, next_time, _vars, &_cur_time, CV_NORMAL);
	if (flag == CV_SUCCESS) return 0;
	else return flag;
}

void ViCaRS::cleanup(void) {
	delete num_send;
	delete num_recv;
	delete send_offset;
	delete recv_offset;
	
	// Free y and abstol vectors
	N_VDestroy_Parallel(_vars);
	N_VDestroy_Parallel(_abs_tol);
	
	// Free integrator memory
	CVodeFree(&_solver);
}

realtype ViCaRS::interaction(BlockGID a, BlockGID b) {
	if (a!=b) return 0.2;
	//if (fabs(a-b)==1) return 0.00001;
	else return 0;
}

void ViCaRS::write_header(FILE *fp) {
	GlobalLocalMap::const_iterator	it;
	BlockGID		gid;
	fprintf(fp, "t ");
	for (it=_global_local_map.begin();it!=_global_local_map.end();++it) {
		gid = it->first;
		fprintf(fp, "x%d v%d h%d F%d ", gid, gid, gid, gid);
	}
	fprintf(fp, "\n");
}

void ViCaRS::write_cur_data(FILE *fp) {
	realtype			xi, vi, hi;
	GlobalLocalMap::const_iterator	it;
	
	fprintf(fp, "%0.7e ", _cur_time);
	for (it=_global_local_map.begin();it!=_global_local_map.end();++it) {
		xi = X(it->first);
		vi = V(it->first);
		hi = H(it->first);
		fprintf(fp, "%14.6e %14.6e %14.6e %14.6e ", xi, vi, hi, F(it->first, vi, hi));
	}
	fprintf(fp, "\n");
}

void ViCaRS::print_stats(void) {
	long int nst, nfe, nni, ncfn, netf, njtv;
	int flag;
	
	flag = CVodeGetNumSteps(_solver, &nst);
	flag = CVodeGetNumRhsEvals(_solver, &nfe);
	flag = CVodeGetNumErrTestFails(_solver, &netf);
	flag = CVodeGetNumNonlinSolvIters(_solver, &nni);
	flag = CVodeGetNumNonlinSolvConvFails(_solver, &ncfn);
	flag = CVSpilsGetNumJtimesEvals(_solver, &njtv);
	
	printf("\nFinal Statistics:\n");
	printf("NumSteps = %-6ld\n", nst);
	printf("NumRhsEvals  = %-6ld\n", nfe);
	printf("NumErrTestFails = %-6ld\n", netf);
	printf("NumNonlinSolvIters = %-6ld\n", nni);
	printf("NumNonlinSolvConvFails = %-6ld\n", ncfn);
	printf("NumJtimesEvals = %-6ld\n", njtv);
}

/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

//#define SLOWNESS_LAW

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
	for (it=sim->begin();it!=sim->end();++it) {
		if (Vth(y,it->second) <= 0 || Hth(y,it->second) <= 0) {
			local_fail = 1;
			break;
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
#ifdef SLOWNESS_LAW
		Hth(ydot,lid) = RCONST(1) - h*v;
#else
		Hth(ydot,lid) = -h*v*log(h*v);
#endif
	}
	
	delete global_x;
	
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
						+ -(sim->param_k(gid)*sim->param_a(gid)/(sim->param_r(gid)*Vth(y,lid))) * Vth(v,lid)
						+ -(sim->param_k(gid)*sim->param_b(gid)/(sim->param_r(gid)*Hth(y,lid))) * Hth(v,lid);
#ifdef SLOWNESS_LAW
		Hth(Jv,lid) =	-Hth(y,lid) * Vth(v,lid)
						+ -Vth(y,lid) * Hth(v,lid);
#else
		Hth(Jv,lid) =	-Hth(y,lid)*(log(Hth(y,lid)*Vth(y,lid)) + 1) * Vth(v,lid)
						+ -Vth(y,lid)*(log(Vth(y,lid)*Hth(y,lid)) + 1) * Hth(v,lid);
#endif
		//std::cerr << lid << " " << gid << " " << Xth(Jv,lid) << " " << Vth(Jv,lid) << " " << Hth(Jv,lid) << std::endl;
	}
	
	return 0;
}