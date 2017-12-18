#include "instance.hpp"
#include "misc.hpp"

// TODO: SMH: does this need to be in here?
extern runType runParam;

instance::instance () {}

void instance::initialize(PowSys *powSys, scenarios *stoc, string inputDir) {
	this->powSys = powSys;
	this->stoc	 = stoc;
	
	solution.allocateMem(powSys->numGen, (int)round(runParam.horizon/runParam.ED_resolution));
	
	DA_load.read(inputDir + "Load/DA.csv");
	ST_load.read(inputDir + "Load/ST.csv");
	RT_load.read(inputDir + "Load/RT.csv");
}
