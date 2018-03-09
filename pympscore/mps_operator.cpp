#include "mps_operator.h"
#include "mps_neighbor.h"
#include "mps_weight.h"

#include <iostream>

#include <Eigen/sparse>

MPS_NAMESPACE_BEGIN;

void calcNumberDensity(
	const WeightFunction & wfun,
	double re,
	const NeighborTable & table,
	EigenRef<VectorXd> values)
{
	//const double re = table.getCutoff();
	assert(table.getCutoff() >= re);

	values.setZero();

	const size_t np = values.size();
	for (int i = 0; i < np; i++) {
		double val = 0;
		for (const auto &neigh : table[i]) {
			if (neigh.distance < re) {
				val += wfun(neigh.distance, re);
			}
		}
		values(i) = val;
	}

}

void tagBoundaryFlag(
	double ntol, 
	EigenRef<const VectorXd>& ndens, 
	EigenRef<VectorXi>& flag)
{
	const size_t npart = ndens.size();
	for (int ipart = 0; ipart < npart; ipart++) {
		if (ndens(ipart) < ntol) {
			flag(ipart) = 1;
		}
		else {
			flag(ipart) = 0;
		}
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
		double pi = values(ipart);

		const bool use_pmin = true;
		if (use_pmin) {
			for (const auto &neigh : table[ipart]) {
				if (neigh.distance < re) {
					pi = std::min(pi, values(neigh.index));
				}
			}
		}


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
				if (use_pmin) {
					pij = values(jpart) - pi;
				}

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

void calcLaplacian(
	const WeightFunction & wfun, 
	double re, double coef, 
	const NeighborTable & table, 
	EigenRef<const VectorXd> phi, 
	EigenRef<VectorXd> lap)
{
	const size_t npart = phi.size();
	for (int ipart = 0; ipart < npart; ipart++) {
		const double iphi = phi(ipart);

		double val = 0;

		for (const auto &neigh : table[ipart]) {
			const int jpart = neigh.index;
			const double rij = neigh.distance;
			if (rij < re) {
				double pij = phi(jpart) - iphi;
				double wij = wfun(rij, re);
				val += pij * wij;
			}
		}

		val *= coef;

		lap(ipart) = val;
	}
}

void solvePres(
	const WeightFunction &wfun,
	double re, double coef,
	const NeighborTable &table,
	EigenRef<const VectorXi> flag,
	EigenRef<VectorXd> rhs,
	EigenRef<VectorXd> sol)
{
	// TODO
	using triplet = Eigen::Triplet<double>;
	using spmat = Eigen::SparseMatrix<double>;



	// build matrix
	const size_t npart = flag.size();

	std::vector<triplet> matbuf;
	matbuf.reserve(npart * 16);

	for (int ipart = 0; ipart < npart; ipart++) {
		if (flag(ipart) == 1) { // surface
			matbuf.emplace_back(ipart, ipart, 1.0);
			rhs(ipart) = 0.0;
		}
		else { // interior
			for (const auto &neigh : table[ipart]) {
				const int jpart = neigh.index;
				const double rij = neigh.distance;
				if (rij < re) {
					double wij = wfun(rij, re);
					double aij = wij * coef;

					matbuf.emplace_back(ipart, ipart, aij);
					//matbuf.emplace_back(ipart, jpart, -aij);
					if (flag(jpart) == 1) {
						// rhs(ipart) += aij * 0;
					}
					else {
						matbuf.emplace_back(ipart, jpart, -aij);
					}
				}
			}
		}
	}

	// build matrix
	spmat mat(npart, npart);
	mat.setFromTriplets(matbuf.begin(), matbuf.end());
	
	//
	Eigen::ConjugateGradient<
		spmat, Eigen::Lower | Eigen::Upper, 
		Eigen::IncompleteCholesky<double> > solver;
	//Eigen::BiCGSTAB<spmat> solver;
	solver.setTolerance(1.0e-12);
	solver.compute(mat);
	sol = solver.solve(rhs);

	std::cout << __FUNCTION__ 
		<< ": cycle=" << solver.iterations() 
		<< ", error=" << solver.error() 
		<< std::endl;

}



MPS_NAMESPACE_END;


////////////////////////////////////////////////////////////////////////////////


#include "pymps_bind.h"


void def_pymod_operator(py::module &mod) {
	using namespace mps;

	mod.def("calcNumDens", &calcNumberDensity, "Number density");
	mod.def("calcSymGrad", &calcSymGradient, "Symmetric gradient pi+pj");
	mod.def("calcLaplace", &calcLaplacian, "Laplacian");
	mod.def("calcPresEos", &calcPresEos, "Pressure by EOS");

	mod.def("tagSurfFlag", &tagBoundaryFlag, "Tag surface particles");
	mod.def("solvePres", &solvePres, "Solve PPE");
}

//--------------------------------------------------------------------------------



