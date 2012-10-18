#include "Equations.h"
#include "RateState.h"

/*
 * Functions called by the solver for advanced RS
 */

bool OrigEqns::values_valid(ViCaRS *sim, N_Vector y) {
	int			local_fail;
    realtype    v, h;
	BlockMap::const_iterator	it;
	
	if (_use_log_spline) return true;
	
	// Check if any velocity or theta values are below 0
	local_fail = 0;
	for (it=sim->begin();it!=sim->end();++it) {
        h = Hth(y,it->first);
        v = Vth(y,it->first);
		if (v <= 0 || h <= 0 || isnan(v) || isnan(h)) {
			local_fail = 1;
			break;
		}
	}
	
	return (!local_fail);
}

bool SimpleEqns::values_valid(ViCaRS *sim, N_Vector y) {
	int			local_fail;
    realtype    h;
	BlockMap::const_iterator	it;
	
	// Check if any velocity or theta values are below 0
	local_fail = 0;
	for (it=sim->begin();it!=sim->end();++it) {
        h = Hth(y,it->first);
		if (h <= 0 || isnan(h)) {
			local_fail = 1;
			break;
		}
	}
	
	return (!local_fail);
}

bool TullisEqns::values_valid(ViCaRS *sim, N_Vector y) {
	return true;
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
		case 3: return "S";
		case 4: return "P";
		default: throw std::exception();
	}
}

std::string TullisEqns::var_name(unsigned int var_num) const {
	switch (var_num) {
		case 0: return "X";
		case 1: return "V";
		case 2: return "Tau";
		default: throw std::exception();
	}
}

realtype OrigEqns::var_value(ViCaRS *sim, unsigned int var_num, BlockGID gid, N_Vector y) {
	switch (var_num) {
		case 0: return Xth(y, gid);
		case 1: return Vth(y, gid);
		case 2: return Hth(y, gid);
		case 3: return F(sim->param_a(gid), sim->param_a(gid), Vth(y, gid), Hth(y, gid));
		/*case 0: return Xth(y, gid)*sim->params().L;
		case 1: return Vth(y, gid)*sim->params().v_ss;
		case 2: return Hth(y, gid)*sim->params().L/sim->params().v_ss;
		case 3: return F(sim->param_a(gid), sim->param_a(gid), Vth(y, gid), Hth(y, gid));*/
		default: throw std::exception();
	}
}

realtype TullisEqns::var_value(ViCaRS *sim, unsigned int var_num, BlockGID gid, N_Vector y) {
	switch (var_num) {
		case 0: return Xth(y, gid);
		case 1: return Vth(y, gid);
		case 2: return Tauth(y, gid);
		default: throw std::exception();
	}
}

realtype SimpleEqns::var_value(ViCaRS *sim, unsigned int var_num, BlockGID gid, N_Vector y) {
	switch (var_num) {
		case 0:
		case 1: return NV_DATA_S(y)[gid*num_equations()+var_num];
		case 2: return _vel[gid];
		case 3: return _elem_stress[gid];
		case 4: return _phase[gid];
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

void TullisEqns::init_block(BlockGID gid, const BlockData &block, N_Vector vars, N_Vector tols) {
	Xth(vars, gid) = block._init_x;
	Vth(vars, gid) = block._init_v;
	Tauth(vars, gid) = block._init_h;
	
	Xth(tols, gid) = block._tol_x;
	Vth(tols, gid) = block._tol_v;
	Tauth(tols, gid) = block._tol_h;
}

void SimpleEqns::init_block(BlockGID gid, const BlockData &block, N_Vector vars, N_Vector tols) {
	_start_time[gid] = 0;
	
	Xth(vars, gid) = block._init_x;
	Hth(vars, gid) = block._init_h;
	
	Xth(tols, gid) = block._tol_x;
	Hth(tols, gid) = block._tol_h;
}

int OrigEqns::init(ViCaRS *sim) {
	return 0;
}

int TullisEqns::init(ViCaRS *sim) {
	return 0;
}

int SimpleEqns::init(ViCaRS *sim) {
	BlockMap::const_iterator	it,jt;
	realtype					sum_load, delta_tau;
	
	_theta_star = sim->params().D_c/sim->params().V_star;
	
	for (it=sim->begin();it!=sim->end();++it) {
		_phase[it->first] = 0;				// Start all blocks in phase 0
		_ss_stress[it->first] = sim->params().sigma*(sim->params().mu_0+(sim->params().A-sim->params().B)*log(sim->params().v_ss/sim->params().V_star));
        std::cerr << _ss_stress[it->first] << std::endl;
		delta_tau = _ss_stress[it->first];
		_v_eq[it->first] = 2*sim->params().beta*delta_tau/sim->params().G;
		_base_stress[it->first] = _elem_stress[it->first] = 0;
		_start_time[it->first] = 0;
	}
	
	for (it=sim->begin();it!=sim->end();++it) {
		sum_load = 0;
		for (jt=sim->begin();jt!=sim->end();++jt) {
			sum_load += sim->interaction(it->first, jt->first);
		}
		sum_load += sim->params().G/sim->params().W;
		_stress_loading[it->first] = sim->params().v_ss*sum_load;
        std::cerr << sum_load << " " << _stress_loading[it->first] << std::endl;
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

// Function to determine when the simulation moves into rupture or long term mode
// TODO: parallelize this function
int TullisEqns::check_for_mode_change(ViCaRS *sim, realtype t, N_Vector y, realtype *gout) {
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
// Term		Phase 0		Phase 1			Phase 2
// dX/dX =	0			0				0
// dX/dH =	0							0
// dH/dX =
// dH/dH =	0
int SimpleEqns::solve_odes(ViCaRS *sim, realtype t, N_Vector y, N_Vector ydot) {
	BlockMap::const_iterator	it, j;
	int             phase_num;
	BlockGID        gid;
	realtype        x, h, v_0=0.001, tau_dot, H, K, q;
	
	// Calculate velocities and stresses on blocks
	for (it=sim->begin();it!=sim->end();++it) {
		gid = it->first;
		switch (_phase[gid]) {
			case 0:
				_vel[gid] = 0;
				_elem_stress[gid] = _base_stress[gid] + (t-_start_time[gid])*_stress_loading[gid];
				break;
			case 1:
				K = -sim->interaction(gid, gid)+1;	// K_T = 1 arbitrarily, need to change this
				H = sim->params().B/sim->params().D_c-K/sim->params().sigma;
				tau_dot = _ss_stress[gid];
				q = 1.0/((1.0/v_0+H*sim->params().sigma/tau_dot)*exp(-tau_dot*(t-_start_time[gid])/(sim->params().A*sim->params().sigma)) - H*sim->params().sigma/tau_dot);
				_vel[gid] = q;
				_elem_stress[gid] = sim->params().sigma*(sim->params().mu_0+sim->params().A*log(_vel[gid]/sim->params().V_star)+sim->params().B*log(Hth(y,gid)/_theta_star));
				break;
			case 2:
				_vel[gid] = _v_eq[gid];
				_elem_stress[gid] = sim->params().sigma*(sim->params().mu_0+sim->params().A*log(_vel[gid]/sim->params().V_star)+sim->params().B*log(Hth(y,gid)/_theta_star));
				break;
		}
	}
	
    std::cerr << gid << " time: " << t/(365.25*86400) << " phase: " << _phase[gid] << " H: " << Hth(y,gid) << " vel: " << _vel[gid] << " stress: " << _elem_stress[gid] << std::endl;
    
	for (it=sim->begin();it!=sim->end();++it) {
		gid = it->first;
		
		phase_num = _phase[gid];
		
		x = Xth(y,gid);
		h = Hth(y,gid);
		
		switch (phase_num) {
			case 0:
				Hth(ydot,gid) = 1;
				break;
			case 1:
			case 2:
				Hth(ydot,gid) = 1.0-h*_vel[gid]/sim->params().D_c;
				break;
			default:
				return -1;
		}
		
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
		
		switch (_phase[gid]) {
			case 0:
				max_val = fmax(max_val, _elem_stress[gid] - _ss_stress[gid]);
				break;
			case 1:
				max_val = fmax(max_val, _vel[gid] - _v_eq[gid]);
				break;
			case 2:
				max_val = fmax(max_val, _ss_stress[gid] - _elem_stress[gid]);
				break;
			default:
				return 1;
		}
	}
	
	gout[0] = max_val;
	
	return 0;
}

void OrigEqns::handle_mode_change(ViCaRS *sim, realtype t, N_Vector y) {
}

void TullisEqns::handle_mode_change(ViCaRS *sim, realtype t, N_Vector y) {
}

void SimpleEqns::handle_mode_change(ViCaRS *sim, realtype t, N_Vector y) {
	BlockMap::const_iterator it;
	BlockGID gid;
    bool        all_phase_0 = true;
	
	for (it=sim->begin();it!=sim->end();++it) {
		gid = it->first;
		if (_phase[gid] == 0 && _elem_stress[gid] >= _ss_stress[gid]) {
			_phase[gid] = 1;
			_start_time[gid] = t;
		} else if (_phase[gid] == 1 && _vel[gid] >= _v_eq[gid]) {
			_phase[gid] = 2;
			_start_time[gid] = t;
		} else if (_phase[gid] == 2 && _ss_stress[gid] >= _elem_stress[gid]) {
			_phase[gid] = 0;
			_start_time[gid] = t;
		}
        if (_phase[gid] != 0) all_phase_0 = false;
	}
}


int TullisEqns::solve_odes(ViCaRS *sim, realtype t, N_Vector y, N_Vector ydot) {
	realtype		friction_force, spring_force, x, v, tau, interact, trm, exptrm, stauss;
    realtype        arg, vs, vpl, scab, ca, cb, coefv, denom, coeft1, coeft2, rigid, strs, xldis;
	unsigned int	i_gid, j_gid;
	BlockMap::const_iterator	it,jt;
	
	// Share X values over all processors to enable correct interaction between blocks
	// TODO: only share X values with processors that need them
	
    vs = 3.0e3; // shear wave speed
    
    vpl = sim->params().v_ss;
    ca = sim->params().A;
    cb = sim->params().B;
    xldis = sim->params().D_c;
    scab = ca - cb;
    
	for (it=sim->begin();it!=sim->end();++it) {
		i_gid = it->first;
		x = Xth(y,i_gid);
		v = Vth(y,i_gid);
		tau = Tauth(y,i_gid);
		
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
		
        strs = 0;
        rigid = 3.0e4;
        denom = 2 * vs * ca + rigid * v;        // (m/s)*(N/m^2)+(N/m^2*m/s)=N/(m*s)
        coefv = 2 * vs * v / denom;             // (m/s)*(m/s)/(N/(m*s))=(m^3/s^3)/N
        coeft1 = 2 * vs * ca / denom;           // (m/s)*(N/m^2)/(N/(m*s))=unitless
        coeft2 = rigid * v / denom;             // (N/m^2)*(m/s)/(N/(m*s))=unitless
        
		Xth(ydot,i_gid) = v;
        arg = vpl/v;
        stauss = -scab * log(arg);              // N/m^2
        exptrm = exp(-(tau-stauss)/cb);         // unitless
        trm = cb * v/xldis * (1-exptrm);        // (N/m^2)*(m/s)/
		Vth(ydot,i_gid) = coefv * (strs + trm); // (m^3/s^3)/N*(N/m^2)=(m/s^3)
        Tauth(ydot,i_gid) = coeft1*strs + coeft2*trm;
	}
	
	return 0;
}

