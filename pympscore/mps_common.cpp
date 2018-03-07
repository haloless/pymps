
#include "mps_common.h"


////////////////////////////////////////////////////////////////////////////////
MPS_NAMESPACE_BEGIN;

static int SpaceDim = 2;

int GetSpaceDim()
{
	return SpaceDim;
}

void SetSpaceDim(int ndim)
{
	SpaceDim = ndim;
}

VectorXd MakeSpaceDimVec()
{
	return VectorXd(SpaceDim);
}

MatrixXd MakeSpaceDimMat()
{
	return MatrixXd(SpaceDim, SpaceDim);
}

MPS_NAMESPACE_END;
////////////////////////////////////////////////////////////////////////////////


#include "pymps_bind.h"

//--------------------------------------------------------------------------------
void def_pymod_common(py::module &mod) 
{
	using namespace mps;
	mod.def("GetSpaceDim", &GetSpaceDim, "Get space dimension");
	mod.def("SetSpaceDim", &SetSpaceDim, "Set space dimension");
}







////////////////////////////////////////////////////////////////////////////////

