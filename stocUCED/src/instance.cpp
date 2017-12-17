#include "instance.hpp"
#include "misc.hpp"

// TODO: SMH: does this need to be in here?
extern runType runParam;

instance::instance () {}

void instance::initialize(PowSys *powSys, scenarios *stoc) {
	this->powSys = powSys;
	this->stoc	 = stoc;
	
	solution.allocateMem(powSys->numGen, 14);
	cout << "Semih: Fix solution mem alloc" << endl;
}

/*
bool instance::getGenStatus(int gen, int period)
{
	if (period < 0)	return true;	// all generators are assumed to be online (for a long time) at t=0.

	if (period >= gen_states[gen].size()) {
		cout << "Generator status has not been determined yet!" << endl;
		return false;
	}

	bool status = true;
	for (int t=0; t<gen_states[gen].size(); t++) {
		if ( (t+1) > period ) {
			status = gen_states[gen][t];
			break;
		}
	}

	return status;
}

void instance::updateGenStatus(ProblemType prob_type, vector< vector<double> > &x)
{
	int nb_subhours = 60 / runParam.ST_resolution;

	for (int g = 0; g < numGen; g++)
	{
		if (prob_type == DayAhead && dayahead_gen[g])	// if the problem was DA-UC, and this generator is a base-load gen, ...
		{
			for (int h = 0; h< runParam.DA_numPeriods; h++) {
				for (int t=0; t<nb_subhours; t++) {
					gen_states[g].push_back( round(x[g][h]) );
				}
			}
		}
		else if (prob_type == ShortTerm && !dayahead_gen[g])
		{
			for (int t = 0; t < runParam.ST_numPeriods; t++) {
				gen_states[g].push_back( round(x[g][t]) );
			}
		}
	}
}*/
