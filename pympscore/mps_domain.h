#pragma once

#include "mps_common.h"

#include <array>
#include <vector>
#include <memory>

////////////////////////////////////////////////////////////////////////////////
MPS_NAMESPACE_BEGIN;

// fwd decl
class BucketCache;
class BucketIter;
class Domain;


//
// TODO periodic???
//
class BucketCache
{
protected:

	VectorXd m_vmin, m_vmax;

	std::vector<double> m_spacing;
	std::vector<int> m_cellnum;

	typedef std::vector<int> bucket_type;
	std::vector<bucket_type> m_buckets;

	MatrixXi m_indices;

public:

	BucketCache(VectorXd const &vmin, VectorXd const &vmax);

	BucketCache(const Domain &domain);

	~BucketCache();

	void defineByRefLength(double reflen);

	void clearAll();

	//void clearAt(int i, int j, int k);

	auto& operator()(int i, int j, int k) { return m_buckets[indexOf(i, j, k)]; }

	auto& operator()(int i, int j) { return m_buckets[indexOf(i, j)]; }

	auto& operator[](const VectorXi &ind) { return m_buckets[indexOf(ind)]; }
	auto& operator[](const Vector3i &ind) { return m_buckets[indexOf(ind)]; }

	const VectorXd& vmin() const { return m_vmin; }
	const VectorXd& vmax() const { return m_vmax; }

	auto& spacing() const { return m_spacing; }
	auto& cellnum() const { return m_cellnum; }


	//----------------------------------------

	void getExtent(Vector3i &lo, Vector3i &hi);

	VectorXi lookup(const VectorXd &pos) const;

	void cachePositions(EigenRef<const MatrixXd> pos);

	VectorXi getCachedIndex(int ipart) const {
		VectorXi ind = m_indices.row(ipart);
		return ind;
	}

	auto& toList() { return m_buckets; }

	std::string toString() const;

protected:
	int indexOf(const VectorXi &ind);
	int indexOf(const Vector3i &ind) { return indexOf(ind(0), ind(1), ind(2)); }
	int indexOf(int i, int j);
	int indexOf(int i, int j, int k);

};

//--------------------------------------------------------------------------------

//
// TODO periodic???
//
class BucketIter
{
protected:

	BucketCache * m_cache = nullptr;

	Vector3i m_range_lo;
	Vector3i m_range_hi;

	int m_nband = 1; // iteration width

	Vector3i m_cen; // center position
	Vector3i m_pos; // current position
	//bool m_ok = false;

public:

	BucketIter(BucketCache &cache, const VectorXi &center, int nband = 1);

	~BucketIter() {}

	// iteration band width
	int getIterWidth() const { return m_nband; }

	//
	const Vector3i& getIterCenter() const { return m_cen; }

	void next();
	bool ok() const;
	void reset();

	const Vector3i& getIterPos() const { return m_pos; }
	const Vector3i& operator()() const { return m_pos; }

	BucketIter& operator++() { next(); return *this; }

protected:

	void setIterWidth(int nband) { m_nband = nband; }

	void setIterCenter(const VectorXi &ind);

	void setValidRange();

};




//--------------------------------------------------------------------------------

class Domain
{
protected:

	// lower and upper corners
	VectorXd vmin, vmax;

	std::unique_ptr<BucketCache> m_cache;


public:

	Domain();

	~Domain();

	////////////////////////////////////////

	const VectorXd& getLowerEnd() const { return vmin; }
	const VectorXd& getUpperEnd() const { return vmax; }
	void setLowerEnd(VectorXd v) { vmin = v; }
	void setUpperEnd(VectorXd v) { vmax = v; }

	// 
	void createBucketCache(double length);

	void destroyCache();

	auto& getCache() { return *m_cache; }

	std::string toString() const;

};





MPS_NAMESPACE_END;
////////////////////////////////////////////////////////////////////////////////



