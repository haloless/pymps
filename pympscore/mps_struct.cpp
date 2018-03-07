#include "mps_struct.h"

#include <iostream>

////////////////////////////////////////////////////////////////////////////////

MPS_NAMESPACE_BEGIN;



FluidData::~FluidData()
{
	//std::cout << __FUNCTION__ << std::endl;
}

void FluidData::alloc() {
	// 
	numCapacity = numParticle;

	alloc_vector_field(position);
	alloc_vector_field(position_prev);

	alloc_vector_field(velocity);
	alloc_vector_field(velocity_prev);

	alloc_scalar_field(pressure);
	alloc_scalar_field(pressure_prev);

	alloc_scalar_field(particleType);
	alloc_scalar_field(particleNumDens);

	alloc_scalar_field(particleFlagBC);
}

void FluidData::alloc_vector_field(MatrixXd & f) {
	f.resize(numCapacity, numDim);
	f.setZero();
}

void FluidData::alloc_scalar_field(VectorXd & f) {
	f.resize(numCapacity);
	f.setZero();
}

void FluidData::alloc_scalar_field(VectorXi & f) {
	f.resize(numCapacity);
	f.setZero();
}

MPS_NAMESPACE_END;
////////////////////////////////////////////////////////////////////////////////

#include "pymps_bind.h"



//--------------------------------------------------------------------------------

void def_pymod_consts(py::module &mod) {
	using namespace mps;

	////////////////////////////////////////
#if 1
	py::enum_<ParticleType>(mod, "ParticleType")
		.value("GHOST", ParticleType::GHOST, "Ghost type")
		//.export_values()
		;
#else
#endif

	////////////////////////////////////////

}

//--------------------------------------------------------------------------------

void def_pymod_struct(py::module &mod) {
	using namespace mps;

	////////////////////////////////////////


	py::class_<FluidData>(mod, "FluidData")
		// init
		.def(py::init<>(), "Init empty")
		.def(py::init<int>(), "Init by dimension")
		// properties
		.def_readwrite("numDim", &FluidData::numDim, "Number of dimensions")
		.def_readwrite("numParticle", &FluidData::numParticle, "Number of particles")
		//
		.def("alloc", &FluidData::alloc, "Alloc fields")

		// fluid fields
		.def("position", &FluidData::getPosition, pyrvp::reference_internal)
		.def("velocity", &FluidData::getVelocity, pyrvp::reference_internal)
		.def("pressure", &FluidData::getPressure, pyrvp::reference_internal)
		.def("numdens", &FluidData::getParticleNumDens, pyrvp::reference_internal)
		;

	////////////////////////////////////////

	py::class_<Timer>(mod, "Timer")
		.def(py::init<>())
		.def_readwrite("dt", &Timer::dt, "Time step")
		.def_readwrite("dt_initial", &Timer::dt_initial)
		.def_readwrite("step", &Timer::step)
		;

	////////////////////////////////////////

	py::class_<FluidProperty>(mod, "FluidProperty")
		.def(py::init<>())
		.def("gravity", &FluidProperty::getGravity, pyrvp::reference_internal)
		;

	////////////////////////////////////////

	py::class_<FluidParameter>(mod, "FluidParameter")
		.def(py::init<>())
		.def_readwrite("nZero", &FluidParameter::nZero, "n0, standard mps particle number density")
		;

	////////////////////////////////////////


}




