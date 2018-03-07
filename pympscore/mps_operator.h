#pragma once

#include "mps_common.h"

MPS_NAMESPACE_BEGIN;

// fwd
class WeightFunction;
class NeighborTable;


// 
void calcNumberDensity(
	const WeightFunction &wfun,
	const NeighborTable &table,
	EigenRef<VectorXd> values);


// 
void calcSymGradient(
	const WeightFunction &wfun,
	double re, double nzero,
	const NeighborTable &table,
	EigenRef<const MatrixXd> coords,
	EigenRef<const VectorXd> values,
	EigenRef<MatrixXd> grad);

//
void calcPresEos(
	EigenRef<const VectorXd> numdens,
	EigenRef<VectorXd> pres,
	double pref, double nzero);



MPS_NAMESPACE_END;







