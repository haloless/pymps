#pragma once

#include "mps_common.h"

MPS_NAMESPACE_BEGIN;

// fwd
class NeighborTable;


class WeightFunction
{
public:

	// evalulate
	virtual double eval(double dist, double cutoff) const = 0;

	// batch evaluate
	virtual void eval_batch(int count, double *weight, const double *dist, double cutoff) const {
		for (int i = 0; i < count; i++) {
			weight[i] = eval(dist[i], cutoff);
		}
	}

	// 
	inline double operator()(double dist, double cutoff) const {
		return eval(dist, cutoff);
	}

};

class StandardWeightFunction : public WeightFunction 
{
public:
	// WeightFunction ‚ð‰î‚µ‚ÄŒp³‚³‚ê‚Ü‚µ‚½
	virtual double eval(double dist, double cutoff) const override {
		if (dist < cutoff) {
			return cutoff / dist - 1;
		}
		else {
			return 0;
		}
	}

	static const std::string& GetDoc();

};

class ModifiedWeightFunction : public WeightFunction 
{
public:
	// WeightFunction ‚ð‰î‚µ‚ÄŒp³‚³‚ê‚Ü‚µ‚½
	virtual double eval(double dist, double cutoff) const override {
		if (dist < cutoff) {
			double q = dist / cutoff;
			return q + 1.0 / q - 2;
		}
		else {
			return 0;
		}
	}

	static const std::string& GetDoc();

};

class ModifiedGradientWeightFunction : public WeightFunction
{
public:
	// WeightFunction ‚ð‰î‚µ‚ÄŒp³‚³‚ê‚Ü‚µ‚½
	virtual double eval(double dist, double cutoff) const override {
		if (dist < cutoff) {
			double q = dist / cutoff;
			return 1.0 / q - q;
		}
		else {
			return 0;
		}
	}
};





MPS_NAMESPACE_END;







