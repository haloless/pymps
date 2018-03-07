
#include "mps_neighbor.h"
#include "mps_domain.h"

#include <iostream>
#include <algorithm>

MPS_NAMESPACE_BEGIN;

NeighborEntry& NeighborList::addNewEntry() {
	emplace_back();
	return back();
}

void NeighborList::sortByDistance() {
	std::sort(begin(), end());
}

////////////////////////////////////////////////////////////////////////////////
NeighborTable::NeighborTable()
{
}

NeighborTable::~NeighborTable()
{
	//std::cout << __FUNCTION__ << std::endl;
}

void NeighborTable::resize(size_t npart)
{
	m_table.resize(npart);
}

void NeighborTable::reserveAll(int cap)
{
	for (auto &x : m_table) {
		x.reserve(cap);
	}
}

void NeighborTable::clearAll() {
	for (auto &x : m_table) {
		x.clear();
	}
}

void NeighborTable::buildFromCache(
	BucketCache & cache, 
	EigenRef<const MatrixXd> coords)
{

	const int ndim = mps::GetSpaceDim();
	const size_t npart = m_table.size();

	// cache range
	const int ncx = cache.cellnum()[0];
	const int ncy = cache.cellnum()[1];
	const int ncz = ndim > 2 ? cache.cellnum()[2] : 1;

	const int nband = 1; // search range

	resize(npart);
	clearAll();

	for (int ipart = 0; ipart < npart; ipart++) {
		VectorXi igrid = cache.getCachedIndex(ipart);
		VectorXd ipos = coords.row(ipart);

#if 0
		int ic = igrid(0);
		int jc = igrid(1);
		int kc = ndim > 2 ? igrid(2) : 0;

		for (int k = kc - nband; k <= kc + nband; k++) {
			if (k < 0 || k >= ncz) continue;
			for (int j = jc - nband; j <= jc + nband; j++) {
				if (j < 0 || j >= ncy) continue;
				for (int i = ic - nband; i <= ic + nband; i++) {
					if (i < 0 || i >= ncx) continue;

					auto &b = cache(i, j, k);
					for (int jpart : b) {
						if (jpart != ipart) {
							VectorXd jpos = coords.row(jpart);
							double dist = (jpos - ipos).norm();
							if (dist < m_cutoff) {
								auto &e = m_table[ipart].addNewEntry();
								e.index = jpart;
								e.distance = dist;
							}
						}
					}
				}
			}
		}
#else
		for (BucketIter iter(cache, igrid); iter.ok(); ++iter) {
			for (int jpart : cache[iter()]) {
				if (jpart != ipart) {
					VectorXd jpos = coords.row(jpart);
					double dist = (jpos - ipos).norm();
					if (dist < m_cutoff) {
						auto &e = m_table[ipart].addNewEntry();
						e.index = jpart;
						e.distance = dist;
					}
				}
			}
		}
#endif
	}
}



////////////////////////////////////////////////////////////////////////////////


MPS_NAMESPACE_END;
////////////////////////////////////////////////////////////////////////////////


#include "pymps_bind.h"




void def_pymod_neighbor(py::module &mod) {
	using namespace mps;

	py::class_<NeighborEntry>(mod, "NeighborEntry")
		.def(py::init<>())
		.def_readwrite("index", &NeighborEntry::index)
		.def_readwrite("distance", &NeighborEntry::distance)
		;

	py::class_<NeighborList>(mod, "NeighborList")
		.def(py::init<>())
		.def("asList", &NeighborList::asList, pyrvp::reference_internal)
		.def("addNewEntry", &NeighborList::addNewEntry, pyrvp::reference_internal)
		.def("sortByDistance", &NeighborList::sortByDistance)
		;

	py::class_<NeighborTable>(mod, "NeighborTable")
		.def(py::init<>())
		.def_property("cutoff", &NeighborTable::getCutoff, &NeighborTable::setCutoff)
		.def("resize", &NeighborTable::resize)
		.def("buildFromCache", &NeighborTable::buildFromCache)
		.def("__len__", &NeighborTable::size)
		.def("__getitem__", &NeighborTable::getNeighborList, pyrvp::reference_internal)
		;

}

//--------------------------------------------------------------------------------



