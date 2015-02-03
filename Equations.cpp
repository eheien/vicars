#include "Equations.h"
#include "RateState.h"

/***************************************************
 Dieterich simplified Equations
 ***************************************************/

int SimpleEqns::init(ViCaRS *sim) {
	BlockMap::const_iterator	it,jt;
	realtype					sum_load, delta_tau;
	
	_theta_star = sim->params().D_c/sim->params().V_star;
	
	for (it=sim->begin();it!=sim->end();++it) {
		_phase[it->first] = 0;				// Start all blocks in phase 0
		_ss_stress[it->first] = sim->params().sigma*(sim->params().mu_0+(sim->params().A-sim->params().B)*log(sim->params().v_ss/sim->params().V_star));
		delta_tau = _ss_stress[it->first];
		_v_eq[it->first] = 2*sim->params().beta*delta_tau/sim->params().G;
		_base_stress[it->first] = _elem_stress[it->first] = 0;
		_start_time[it->first] = 0;
        //std::cerr << "Steady state stress: " << _ss_stress[it->first] << " v_eq: " << _v_eq[it->first] << std::endl;
	}
	
	for (it=sim->begin();it!=sim->end();++it) {
		sum_load = 0;
		for (jt=sim->begin();jt!=sim->end();++jt) {
			sum_load += sim->interaction(it->first, jt->first);
		}
		sum_load += sim->params().G/sim->params().W;
		_stress_loading[it->first] = sim->params().v_ss*sum_load;
        //std::cerr << sum_load << " " << _stress_loading[it->first] << std::endl;
	}
	
	return 0;
}

void SimpleEqns::init_block(SimParams params, BlockGID gid, BlockData &block, N_Vector vars, N_Vector tols) {
	_start_time[gid] = 0;
	
	Xth(vars, gid) = 0;
	Hth(vars, gid) = 0.01;
	
	Xth(tols, gid) = 1e-6;
	Hth(tols, gid) = 1e-6;
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
	
    //std::cerr << gid << " time: " << t/(365.25*86400) << " phase: " << _phase[gid] << " H: " << Hth(y,gid) << " vel: " << _vel[gid] << " stress: " << _elem_stress[gid] << std::endl;
    
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

