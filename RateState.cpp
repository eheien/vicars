#include "RateState.h"

ViCaRS::ViCaRS(unsigned int total_num_blocks) : _num_global_blocks(total_num_blocks)
{
    _use_slowness_law = true;
	_G = 3.0e10;
}

int ViCaRS::add_block(const BlockGID &id, const BlockData &block_data) {
	_bdata[id] = block_data;
	
	return 0;
}

int ViCaRS::init(SimEquations *eqns) {
	unsigned int	num_local, i;
	int				flag;
	
	_eqns = eqns;
	_num_equations = eqns->num_equations();
	
	flag = fill_greens_matrix();
	if (flag) return flag;
	
	num_local = (unsigned int)_bdata.size();
	_vars = N_VNew_Serial(_num_global_blocks*_num_equations);
	if (_vars == NULL) return 1;
	
	// Init CVodes structures
	flag = _eqns->init(this);
	if (flag) return flag;
    
	for (i=0;i<_solvers.size();++i) _solvers[i]->init_solver(this);
	
	_cur_solver = 0;
	_cur_time = 0;
	
	return 0;
}

int ViCaRS::advance(void) {
	bool	next_solver;
	int		flag;
	
	flag = _solvers[_cur_solver]->advance(this, _vars, _cur_time, _cur_time, next_solver);
	if (next_solver) {
		// Reinit the next solver
		_cur_solver = (_cur_solver+1)%_solvers.size();
		flag = _solvers[_cur_solver]->reinit_solver(_cur_time, _vars);
		if (flag != 0) return flag;
	}
	
	return flag;
}

void ViCaRS::cleanup(void) {
	unsigned int		i;
	
	// Free y and abstol vectors
	N_VDestroy_Serial(_vars);

	for (i=0;i<_solvers.size();++i) delete _solvers[i];
}

realtype ViCaRS::interaction(BlockGID i, BlockGID j) { return 0; /*greens_matrix->get(i,j);*/ };

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
	fprintf(fp, "%0.2e\t\t%0.2e(%d)\t%0.2e(%d)\n", _cur_time, min_v, min_v_gid, max_v, max_v_gid);
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
			/*if (i==j) { greens_matrix.setVal(i,j,scale*-0.53*G/L); }
			else {
				r = abs(i-j);
				greens_matrix.setVal(i,j,scale*0.98*G/(L*pow(((r/L)-1),3)));
			}*/
		}
	}
	
	return 0;
}

