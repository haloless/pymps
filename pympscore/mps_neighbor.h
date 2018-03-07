#pragma once

#include "mps_common.h"

#include <vector>


MPS_NAMESPACE_BEGIN;

// fwd decl
class Domain;
class BucketCache;



struct NeighborEntry
{
	int index;
	double distance;


	bool operator<(const NeighborEntry &rhs) const {
		return this->distance < rhs.distance;
	}

};

class NeighborList : public std::vector<NeighborEntry>
{
protected:
	typedef NeighborEntry entry_type;
	typedef std::vector<NeighborEntry> list_type;

public:
	NeighborList() : list_type() {}


	list_type& asList() { return *this; }

	NeighborEntry& addNewEntry();

	void sortByDistance();
};


class NeighborTable
{
protected:

	double m_cutoff = 0;

	std::vector<NeighborList> m_table;


public:

	NeighborTable();

	~NeighborTable();

	double getCutoff() const { return m_cutoff; }
	void setCutoff(double cutoff) { m_cutoff = cutoff; }


	const NeighborList & operator[](int i) const { return m_table[i]; }

	NeighborList & getNeighborList(int i) { return m_table[i]; }

	int size() const { return (int)m_table.size(); }

	void resize(size_t npart);

	void reserveAll(int cap);

	void clearAll();

	void buildFromCache(BucketCache &cache, EigenRef<const MatrixXd> coords);
};




MPS_NAMESPACE_END;





