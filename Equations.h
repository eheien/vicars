#include "BlockData.h"
#include "Spline.h"

#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>

#include <string>

#ifndef _VICARS_EQUATIONS_H_
#define _VICARS_EQUATIONS_H_

class ViCaRS;

class SimEquations {
public:
	virtual ~SimEquations(void) {};
	virtual int init(ViCaRS *sim) = 0;
	virtual unsigned int num_equations(void) const = 0;
	virtual unsigned int num_outputs(void) const = 0;
	
	virtual void init_block(BlockGID gid, const BlockData &block, N_Vector vars, N_Vector tols) = 0;
	
	virtual std::string var_name(unsigned int var_num) const = 0;
	virtual realtype var_value(ViCaRS *sim, unsigned int var_num, BlockGID gid, N_Vector y) = 0;
	
	virtual int solve_odes(ViCaRS *sim, realtype t, N_Vector y, N_Vector ydot) = 0;
	virtual bool has_jacobian(void) = 0;
	virtual int jacobian_times_vector(ViCaRS *sim, N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, N_Vector tmp) = 0;
	virtual int check_for_mode_change(ViCaRS *sim, realtype t, N_Vector y, realtype *gout) = 0;
	virtual void handle_mode_change(ViCaRS *sim, realtype t, N_Vector y) = 0;
	virtual bool values_valid(ViCaRS *sim, N_Vector y) = 0;
};

class OrigEqns : public SimEquations {
private:
	LogSpline				log_approx;
    BasicSpline             blog;
	
	// Whether to use the log spline approximation or not
	bool                    _use_log_spline;
	
public:
	OrigEqns(void) : log_approx(-40, 2, 1e12, 1, 15), _use_log_spline(false), blog(1e-5) {};
	virtual ~OrigEqns(void) {};
	
	virtual int init(ViCaRS *sim);
	virtual unsigned int num_equations(void) const { return 3; };
	virtual unsigned int num_outputs(void) const { return 4; };
	
	virtual void init_block(BlockGID gid, const BlockData &block, N_Vector vars, N_Vector tols);
	
	virtual std::string var_name(unsigned int var_num) const;
	virtual realtype var_value(ViCaRS *sim, unsigned int var_num, BlockGID gid, N_Vector y);
	
	virtual int solve_odes(ViCaRS *sim, realtype t, N_Vector y, N_Vector ydot);
	virtual bool has_jacobian(void) { return true; };
	virtual int jacobian_times_vector(ViCaRS *sim, N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, N_Vector tmp);
	virtual int check_for_mode_change(ViCaRS *sim, realtype t, N_Vector y, realtype *gout);
	virtual void handle_mode_change(ViCaRS *sim, realtype t, N_Vector y);
	virtual bool values_valid(ViCaRS *sim, N_Vector y);
    
	realtype &Xth(N_Vector y, BlockGID bnum) { return NV_Ith_S(y,bnum*3+0); };
	realtype &Vth(N_Vector y, BlockGID bnum) { return NV_Ith_S(y,bnum*3+1); };
	realtype &Hth(N_Vector y, BlockGID bnum) { return NV_Ith_S(y,bnum*3+2); };
	
	realtype log_func(realtype x) const {
		if (_use_log_spline) return blog(x);
		else return log(x);
	}
	
	realtype log_deriv(realtype x) const {
		if (_use_log_spline) return blog.dx(x);
		else return 1/x;
	}
	
	realtype F(realtype a, realtype b, realtype v, realtype h) {
		return 1 + (v < 0 ? -1 : 1)*a*log_func(v) + b*log_func(h);
	}
};

class SimpleEqns : public SimEquations {
private:
	// Phase of each element
	std::map<BlockGID, int> _phase;
	
	std::map<BlockGID, realtype>	_ss_stress, _stress_loading, _start_time, _vel, _v_eq;
	std::map<BlockGID, realtype>	_base_stress, _elem_stress;
	
	realtype		_theta_star;
	
public:
	SimpleEqns(void) {};
	virtual ~SimpleEqns(void) {};
	virtual int init(ViCaRS *sim);
	virtual unsigned int num_equations(void) const { return 2; };
	virtual unsigned int num_outputs(void) const { return 5; };
	
	virtual void init_block(BlockGID gid, const BlockData &block, N_Vector vars, N_Vector tols);
	
	virtual std::string var_name(unsigned int var_num) const;
	virtual realtype var_value(ViCaRS *sim, unsigned int var_num, BlockGID gid, N_Vector y);
	
	virtual int solve_odes(ViCaRS *sim, realtype t, N_Vector y, N_Vector ydot);
	virtual bool has_jacobian(void) { return false; };
	virtual int jacobian_times_vector(ViCaRS *sim, N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, N_Vector tmp) { return -1; };
	virtual int check_for_mode_change(ViCaRS *sim, realtype t, N_Vector y, realtype *gout);
	virtual void handle_mode_change(ViCaRS *sim, realtype t, N_Vector y);
	virtual bool values_valid(ViCaRS *sim, N_Vector y);
	
	realtype &Xth(N_Vector y, BlockGID bnum) { return NV_Ith_S(y,bnum*2+0); };
	realtype &Hth(N_Vector y, BlockGID bnum) { return NV_Ith_S(y,bnum*2+1); };
};

class TullisEqns : public SimEquations {
private:
public:
	TullisEqns(void) {};
	virtual ~TullisEqns(void) {};
	virtual int init(ViCaRS *sim);
	virtual unsigned int num_equations(void) const { return 3; };
	virtual unsigned int num_outputs(void) const { return 3; };
	
	virtual void init_block(BlockGID gid, const BlockData &block, N_Vector vars, N_Vector tols);
	
	virtual std::string var_name(unsigned int var_num) const;
	virtual realtype var_value(ViCaRS *sim, unsigned int var_num, BlockGID gid, N_Vector y);
	
	virtual int solve_odes(ViCaRS *sim, realtype t, N_Vector y, N_Vector ydot);
	virtual bool has_jacobian(void) { return false; };
	virtual int jacobian_times_vector(ViCaRS *sim, N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, N_Vector tmp) { return -1; };
	virtual int check_for_mode_change(ViCaRS *sim, realtype t, N_Vector y, realtype *gout);
	virtual void handle_mode_change(ViCaRS *sim, realtype t, N_Vector y);
	virtual bool values_valid(ViCaRS *sim, N_Vector y);
	
	realtype &Xth(N_Vector y, BlockGID bnum) { return NV_Ith_S(y,bnum*3+0); };
	realtype &Vth(N_Vector y, BlockGID bnum) { return NV_Ith_S(y,bnum*3+1); };
	realtype &Tauth(N_Vector y, BlockGID bnum) { return NV_Ith_S(y,bnum*3+2); };
};

#endif
