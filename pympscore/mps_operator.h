#pragma once

#include "mps_common.h"

MPS_NAMESPACE_BEGIN;

// fwd
class WeightFunction;
class NeighborTable;


class OperatorBase {
public:

};


// 
void calcNumberDensity(
	const WeightFunction &wfun,
	double re,
	const NeighborTable &table,
	EigenRef<VectorXd> values);

//
void tagBoundaryFlag(
	double ntol,
	EigenRef<const VectorXd> &ndens,
	EigenRef<VectorXi> &flag);


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

//
void calcLaplacian(
	const WeightFunction &wfun, 
	double re, double coef,
	const NeighborTable &table,
	EigenRef<const VectorXd> phi,
	EigenRef<VectorXd> lap);

//
void solvePres(
	const WeightFunction &wfun,
	double re, double coef,
	const NeighborTable &table, 
	EigenRef<const VectorXi> flag,
	EigenRef<VectorXd> rhs,
	EigenRef<VectorXd> sol);




MPS_NAMESPACE_END;







