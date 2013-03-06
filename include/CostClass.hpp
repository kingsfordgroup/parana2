#ifndef COSTCLASS_HPP
#define COSTCLASS_HPP

#include <unordered_map>
#include <vector>
#include <limits>
#include <boost/multiprecision/gmp.hpp>
#include <boost/pool/pool_alloc.hpp>

#include <cln/rational.h>
#include <cln/integer.h>
#include <cln/float.h>
#include <cln/real.h>

#include "CountedDerivation.hpp"
#include "FHGraph.hpp"
#include "ParanaCommon.hpp"

using std::unordered_map;
using std::vector;
using cln::cl_I;
using boost::multiprecision::mpz_int;
using boost::multiprecision::mpq_rational;

/** Holds the relevant information, for the eager algorithm, about
 * derivations using a particular edge.
 *
 * For each edge and each cost class, there will be an associated
 * EdgeDerivInfoEager instance.  This class will hold the information
 * about the count attributed to each edge, as well as the frontier of
 * the derivations using this edge.
 */
class EdgeDerivInfoEager{

    // typedefs
    // typedef mpz_int BigInt;
    typedef cl_I BigInt;


public:
    typedef std::vector<size_t, CustomAllocator<size_t>> FrontierT;
public:
    EdgeDerivInfoEager( ) : _count(BigInt(0)), _frontier( vector<size_t, CustomAllocator<size_t>>(0,0) ) {}
    EdgeDerivInfoEager( const BigInt& count,
                        const std::vector<size_t, CustomAllocator<size_t>>& frontier ): _count(count),
                        _frontier(std::move(frontier)) {}

    void updateWithDeriv( const CountedDerivation& d ) {
        updateCount( d.count );
        updateFrontier( d.bp );
    }

    void updateFrontier( const std::vector<size_t, CustomAllocator<size_t>>& bp ) {
        assert(bp.size() == _frontier.size());
        for( size_t ti = 0; ti < _frontier.size(); ++ti ) {
            _frontier[ti] = std::max(bp[ti], _frontier[ti]);
        }
    }

    void updateCount( const BigInt& nc ) { _count += nc; }
    const BigInt& count() { return _count; }
    const std::vector<size_t, CustomAllocator<size_t>>& frontier() { return _frontier; }
    void freeDerivations() { _frontier.clear(); }

private:
    FrontierT _frontier;
    BigInt _count;
};

// lazy

class EdgeDerivInfoLazy{
    // typedefs

public:
    //typedef mpz_int BigInt;
    typedef cl_I BigInt;
    typedef std::vector< std::vector<size_t, CustomAllocator<size_t>> > FrontierT;
public:
    EdgeDerivInfoLazy( ) : _count(BigInt(0)), _frontier( vector<vector<size_t, CustomAllocator<size_t>>>(0,vector<size_t, CustomAllocator<size_t>>()) ) {}
    EdgeDerivInfoLazy( const BigInt& count,
                   const std::vector<size_t, CustomAllocator<size_t>>& ibp ): _count(count), _frontier({ibp}) {}

    void updateWithDeriv( const CountedDerivation& d ) {
        updateCount( d.count );
        updateFrontier( d.bp );
    }

    void updateFrontier( const std::vector<size_t, CustomAllocator<size_t>>& bp ) {
        _frontier.push_back(bp);
    }

    void updateCount( const BigInt& nc ) { _count += nc; }

    void freeDerivations() { _frontier.clear(); }

    const BigInt& count() { return _count; }
    const std::vector<std::vector<size_t, CustomAllocator<size_t>>>& frontier() { return _frontier; }

private:
    FrontierT _frontier;
    BigInt _count;
};

template <typename EdgeDerivInfoT>
class CostClass {

    typedef size_t edgeID_t;
    //typedef mpz_int BigInt;
    typedef cl_I BigInt;


public:
    CostClass(): _usedEdges(unordered_map<edgeID_t, EdgeDerivInfoT>()), _cost(-std::numeric_limits<double>::max()) {}
    CostClass(double cost): _usedEdges(unordered_map<edgeID_t, EdgeDerivInfoT>()), _cost(cost) {}

    /** Is the edge eid used by this derivation */
    bool hasEdge( const edgeID_t& eid ) { return _usedEdges.find(eid) != _usedEdges.end(); }

    /**
     * Append this derivation to the cost class, and add it's number
     * of ways of obtaining the given cost solution to the total count
     */
    void appendToCount( const CountedDerivation& deriv ) {
        auto numSln = deriv.count;
        // update the total # of solutions for this cost class
        _total += numSln;
        // the edge used by this derivation
        auto e = deriv.edge;

        // If we have a solution using this edge, then update
        // the highest cost class used by each of the tail vertices
        if ( hasEdge(e) ) {
            // For each tail node, possibly update the largest used
            // cost class
            _usedEdges[e].updateWithDeriv( deriv );
        } else {
            _usedEdges[e] = EdgeDerivInfoT( numSln, deriv.bp  );
        }
    }

    const typename EdgeDerivInfoT::FrontierT& getEdgeFrontier( const edgeID_t& eid ) {
        return _usedEdges[eid].frontier();
    }

    // lazy
    //const vector<vector<size_t>>& getEdgeFrontier( const edgeID_t& eid ) {
    //return _usedEdges[eid].frontier();
    //}

    //eager
    /*
    const vector<size_t>& getEdgeFrontier( const edgeID_t& eid ) {
        return _usedEdges[eid].frontier();
    }
    */
    const vector<edgeID_t> usedEdges() {
       vector<edgeID_t> ue; ue.reserve(_usedEdges.size());
       for( auto it = _usedEdges.begin(); it != _usedEdges.end(); ++it ) {
           if( it->first < std::numeric_limits<size_t>::max() ) {
               ue.push_back(it->first);
           }
       }
       return ue;
    }

    /** Return the probability of edge eid under this score class
     */
    double edgeProb( const edgeID_t& eid ) {
        BigInt edgeCount(0);
        if( hasEdge(eid) ) {
          edgeCount = _usedEdges[eid].count();
        }
        return double_approx(edgeCount / _total);
        /*
        rationalCache_ = edgeCount;
        rationalCache_ /= _total;
        return rationalCache_.convert_to<double>();
        */
    }

   /** Return the total number of derivations of this cost class
   */
   BigInt& total() { return _total; }

   double cost() { return _cost; }

   /** Return the total number of derivations of this cost class
    */
   const BigInt& total() const { return _total; }

   double cost() const { return _cost; }


    void freeDerivations() {
        for ( auto& ueit : _usedEdges ) {
            ueit.second.freeDerivations();
        }
    }
private:
    /** The edges used to obtain this cost class
     * A map from the edge id to the least and greatest cost classes
     * used by each tail node
     */
    unordered_map<edgeID_t, EdgeDerivInfoT> _usedEdges;
    double _cost;
    BigInt _total;
    //mpq_rational rationalCache_;
};

#endif // COSTCLASS_HPP
