#include "mps_weight.h"
#include "mps_neighbor.h"

MPS_NAMESPACE_BEGIN;

const std::string & StandardWeightFunction::GetDoc()
{
	static const std::string doc = 
		"Standard weight function by Koshizuka:\n"
		"w(r,re) = re/r - 1"
		;
	return doc;
}

const std::string & ModifiedWeightFunction::GetDoc()
{
	static const std::string doc =
		"Modified weight function used by e.g. Yamada:\n"
		"w(r,re) = re/r + r/re - 2"
		;
	return doc;
}




////////////////////////////////////////////////////////////////////////////////





////////////////////////////////////////////////////////////////////////////////



MPS_NAMESPACE_END;
////////////////////////////////////////////////////////////////////////////////


#include "pymps_bind.h"


void def_pymod_weight(py::module &mod) {
	using namespace mps;

	auto def_wfun = [&](auto wfun, const char *name) {
		using T = decltype(wfun);
		py::class_<T, WeightFunction>(mod, name)
			.def(py::init<>())
			.def("eval", &T::eval)
			.def("__call__", &T::operator())
			.def("eval_vec", py::vectorize(&T::eval))
			;
	};

	// this is the empty base class
	py::class_<WeightFunction>(mod, "WeightFunction");

	// actual weight functions
	def_wfun(StandardWeightFunction(), "StandardWeightFunction");
	def_wfun(ModifiedWeightFunction(), "ModifiedWeightFunction");
	def_wfun(ModifiedGradientWeightFunction(), "ModifiedGradientWeightFunction");

}

//--------------------------------------------------------------------------------

