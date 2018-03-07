#include "mps_domain.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <iterator>
#include <algorithm>

#include <cmath>
#include <cassert>

////////////////////////////////////////////////////////////////////////////////
MPS_NAMESPACE_BEGIN;


BucketCache::BucketCache(VectorXd const &vmin, VectorXd const &vmax) {
	assert(vmin.size() == vmax.size());

	m_vmin = vmin;
	m_vmax = vmax;

	const size_t ndim = vmin.size();
	assert(ndim == mps::GetSpaceDim());

	m_spacing.resize(ndim);
	m_cellnum.resize(ndim);
}

BucketCache::BucketCache(const Domain & domain)
	: BucketCache(domain.getLowerEnd(), domain.getUpperEnd())
{
}

BucketCache::~BucketCache()
{
	//std::cout << __FUNCTION__ << std::endl;
}

void BucketCache::defineByRefLength(double reflen)
{
	VectorXd vlen = m_vmax - m_vmin;

	const int ndim = mps::GetSpaceDim();
	for (int idir = 0; idir < ndim; idir++) {
		m_cellnum[idir] = (int)(vlen(idir) / reflen);
		m_spacing[idir] = vlen(idir) / m_cellnum[idir];
	}

	int nbucket = 1;
	for (int idir = 0; idir < ndim; idir++) {
		nbucket *= m_cellnum[idir];
	}

	m_buckets.resize(nbucket);

	for (auto &b : m_buckets) {
		b.reserve(24);
	}
}

void BucketCache::clearAll()
{
	for (auto &b : m_buckets) {
		b.clear();
	}
}

void BucketCache::getExtent(Vector3i& lo, Vector3i& hi)
{
	// low end is always zero
	lo.fill(0);

	// copy high end
	hi.fill(0);
	std::transform(m_cellnum.cbegin(), m_cellnum.cend(), hi.data(), 
		[](int ncell) -> int { return ncell - 1; });
}

VectorXi BucketCache::lookup(const VectorXd & pos) const
{
	const int ndim = mps::GetSpaceDim();
	VectorXi ind(ndim);
	
	for (int dir = 0; dir < ndim; dir++) {
		ind(dir) = (int)floor((pos(dir) - m_vmin(dir)) / m_spacing[dir]);
	}

	return ind;
}

void BucketCache::cachePositions(EigenRef<const MatrixXd> coords)
{
	const int ndim = mps::GetSpaceDim();
	assert(coords.cols() == ndim);


	const size_t npart = coords.rows();
	//std::cout << "np=" << npart << std::endl;
	//std::cout << "ndim=" << ndim << std::endl;

	this->clearAll();
	m_indices.resize(npart, ndim);
	//std::cout << m_indices.size() << std::endl;

	for (int i = 0; i < npart; i++) {
		VectorXd pos = coords.row(i);
		VectorXi ind = lookup(pos);

		m_buckets[indexOf(ind)].push_back(i);
		m_indices.row(i) = ind;
	}
}

std::string BucketCache::toString() const
{
	const int ndim = mps::GetSpaceDim();

	std::ostringstream oss;
	oss << "<BucketCache: ";
	oss << "cellnum=(";
	std::copy(m_cellnum.begin(), m_cellnum.end(), std::ostream_iterator<int>(oss, " "));
	oss << "), ";
	oss << "spacing=(";
	std::copy(m_spacing.begin(), m_spacing.end(), std::ostream_iterator<double>(oss, " "));
	oss << "), ";
	oss << ">";

	return oss.str();
}

int BucketCache::indexOf(const VectorXi & ind)
{
	int idx = ind(0) + ind(1)*m_cellnum[0];
	if (mps::GetSpaceDim() > 2) {
		idx += ind(2)*m_cellnum[0] * m_cellnum[1];
	}
	return idx;
}

int BucketCache::indexOf(int i, int j)
{
	return i + j * m_cellnum[0];
}

int BucketCache::indexOf(int i, int j, int k)
{
	int idx = i + j * m_cellnum[0];
	if (mps::GetSpaceDim() > 2) {
		idx += k * m_cellnum[0] * m_cellnum[1];
	}
	return idx;
}

//void BucketCache::clearAt(int i, int j, int k)
//{
//}

////////////////////////////////////////////////////////////////////////////////


BucketIter::BucketIter(BucketCache & cache, const VectorXi & center, int nband)
	: m_cache(&cache), m_nband(nband)
{
	setIterCenter(center);
	setValidRange();

	// ready to use
	reset();
}



void BucketIter::setIterCenter(const VectorXi & ind)
{
	const int ndim = mps::GetSpaceDim();
	m_cen.setZero();
	for (int dir = 0; dir < ndim; dir++) {
		m_cen(dir) = ind(dir);
	}
}

void BucketIter::next()
{
	m_pos(0)++;
	if (m_pos(0) > m_range_hi(0)) {
		m_pos(0) = m_range_lo(0);
		m_pos(1)++;
		if (m_pos(1) > m_range_hi(1)) {
			m_pos(1) = m_range_lo(1);
			m_pos(2)++;
		}
	}
}

bool BucketIter::ok() const
{
	return m_pos(0) <= m_range_hi(0) && m_pos(1) <= m_range_hi(1) && m_pos(2) <= m_range_hi(2);
}

void BucketIter::reset()
{
	m_pos = m_range_lo;
}

void BucketIter::setValidRange() {
	assert(m_cache);

	// cache extent
	Vector3i ext_lo, ext_hi;
	m_cache->getExtent(ext_lo, ext_hi);

	// iter range
	m_range_lo = m_cen - Vector3i::Ones() * m_nband;
	m_range_hi = m_cen + Vector3i::Ones() * m_nband;

	// clip if necessary
	for (int dir = 0; dir < 3; dir++) {
		m_range_lo(dir) = std::max(m_range_lo(dir), ext_lo(dir));
		m_range_hi(dir) = std::min(m_range_hi(dir), ext_hi(dir));
	}
}


////////////////////////////////////////////////////////////////////////////////

Domain::Domain() 
{
	const int ndim = mps::GetSpaceDim();
	vmin.resize(ndim);
	vmax.resize(ndim);
}

Domain::~Domain()
{
	//std::cout << __FUNCTION__ << std::endl;
}

void Domain::createBucketCache(double length)
{
	m_cache.reset(new BucketCache(vmin, vmax));
	m_cache->defineByRefLength(length);
}

void Domain::destroyCache()
{
	m_cache.reset();
}

std::string Domain::toString() const
{
	//static const std::string name = "Domain";
	std::ostringstream oss;
	oss << "<Domain: lo=(" << vmin.transpose() << "), hi=(" << vmax.transpose() << ")>";
	return oss.str();
}

MPS_NAMESPACE_END;
////////////////////////////////////////////////////////////////////////////////


#include "pymps_bind.h"

void def_pymod_domain(py::module &mod) {
	using namespace mps;

	////////////////////////////////////////

	py::class_<BucketCache>(mod, "BucketCache")
		.def(py::init<const Domain&>())
		.def("__repr__", &BucketCache::toString)
		.def("define", &BucketCache::defineByRefLength, py::arg("ref_len"))
		.def_property_readonly("spacing", &BucketCache::spacing)
		.def_property_readonly("cellnum", &BucketCache::cellnum)
		.def("cachePositions", &BucketCache::cachePositions)
		.def("toList", &BucketCache::toList)
		;

	py::class_<Domain>(mod, "Domain")
		.def(py::init<>())
		.def("__repr__", &Domain::toString)
		.def_property("vlo", &Domain::getLowerEnd, &Domain::setLowerEnd)
		.def_property("vhi", &Domain::getUpperEnd, &Domain::setUpperEnd)
		.def("createCache", &Domain::createBucketCache)
		.def("getCache", &Domain::getCache, pyrvp::reference_internal)
		;
}

//--------------------------------------------------------------------------------

