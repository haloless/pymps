#include "mps_operator.h"
#include "mps_neighbor.h"
#include "mps_weight.h"

MPS_NAMESPACE_BEGIN;

void calcNumberDensity(
	const WeightFunction & wfun,
	const NeighborTable & table,
	EigenRef<VectorXd> values)
{
	const double re = table.getCutoff();

	values.setZero();

	const size_t np = values.size();
	for (int i = 0; i < np; i++) {
		double val = 0;
		for (const auto &neigh : table[i]) {
			val += wfun(neigh.distance, re);
		}
		values(i) = val;
	}

}




void calcSymGradient(
	const WeightFunction &wfun,
	double re, double nzero,
	const NeighborTable &table,
	EigenRef<const MatrixXd> coords,
	EigenRef<const VectorXd> values,
	EigenRef<MatrixXd> grad) 
{
	const int ndim = mps::GetSpaceDim();
	assert(coords.cols() == ndim);
	assert(grad.cols() == ndim);

	grad.setZero();

	const size_t npart = coords.rows();
	for (int ipart = 0; ipart < npart; ipart++) {
		const VectorXd ipos = coords.row(ipart);
		const double pi = values(ipart);

		VectorXd igrad(ndim); 
		igrad.setZero();

		for (const auto &neigh : table[ipart]) {
			const int &jpart = neigh.index;
			const double rij = neigh.distance;
			if (rij < re) {
				// rj - ri
				VectorXd xij = coords.row(jpart);
				xij -= ipos;

				double wij = wfun(rij, re);

				double pij = values(jpart) + pi;

				igrad += pij / (rij*rij) * wij * xij;
			}
		}

		grad.row(ipart) = ((double)ndim / nzero) * igrad;
	}
}

void calcPresEos(
	EigenRef<const VectorXd> numdens, 
	EigenRef<VectorXd> pres, 
	double pref, double nzero)
{
	const int gamma = 7;

	pres.fill(0);

	const size_t npart = numdens.size();
	for (int ipart = 0; ipart < npart; ipart++) {
		double nrel = numdens(ipart) / nzero;

		if (nrel >= 1) {
			pres(ipart) = pref * (pow(nrel, gamma) - 1) / gamma;
		}
	}



}

MPS_NAMESPACE_END;


////////////////////////////////////////////////////////////////////////////////


#include "pymps_bind.h"


void def_pymod_operator(py::module &mod) {
	using namespace mps;

	mod.def("calcNumDens", &calcNumberDensity);
	mod.def("calcSymGrad", &calcSymGradient);
	mod.def("calcPresEos", &calcPresEos);
}

//--------------------------------------------------------------------------------



