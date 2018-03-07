#pragma once

#include "mps_common.h"

#include <vector>
#include <memory>

MPS_NAMESPACE_BEGIN;

// fwd
class BucketCache;
class Domain;


enum class ParticleType
{
	GHOST = -1,

};



struct FluidData
{

	int numDim = 0; // dimension
	int numParticle = 0; // particle number
	int numCapacity = 0; // capacity of max number

	// coordinate
	MatrixXd position;
	MatrixXd position_prev;

	// velocity
	MatrixXd velocity;
	MatrixXd velocity_prev;

	// pressure
	VectorXd pressure;
	VectorXd pressure_prev;

	// type
	VectorXi particleType;

	// number density
	VectorXd particleNumDens;

	// 
	VectorXi particleFlagBC;

public:

	FluidData() = default;

	FluidData(int ndim) : numDim(ndim) {}

	~FluidData();

	void alloc();

	//////////

	auto& getPosition() { return position; }
	auto& getPrevPosition() { return position_prev; }

	auto& getVelocity() { return velocity; }
	auto& getPrevVelocity() { return velocity_prev; }

	auto& getPressure() { return pressure; }

	auto& getParticleType() { return particleType; }
	auto& getParticleNumDens() { return particleNumDens; }
	auto& getParticleFlagBC() { return particleFlagBC; }

protected:
	void alloc_vector_field(MatrixXd &f);
	void alloc_scalar_field(VectorXd &f);
	void alloc_scalar_field(VectorXi &f);
};


struct Timer 
{
	double dt;
	double dt_initial;

	int step;


};

struct FluidProperty 
{
	VectorXd gravity;




	auto& getGravity() { return gravity; }
};

struct FluidParameter
{
	double nZero;
};





MPS_NAMESPACE_END;




