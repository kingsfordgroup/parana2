#ifndef COSTCLASS_HPP
#define COSTCLASS_HPP

#include <unordered_map>
#include <vector>
#include <limits>

#include <cln/rational.h>
#include <cln/integer.h>
#include <cln/float.h>
#include <cln/real.h>

#include "CountedDerivation.hpp"
#include "FHGraph.hpp"

using std::unordered_map;
using std::vector;
using cln::cl_I;

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
public:
    typedef std::vector<size_t> FrontierT;
public:
    EdgeDerivInfoEager( ) : _count(cl_I(0)), _frontier( vector<size_t>(0,0) ) {}
    EdgeDerivInfoEager( const cl_I& count,
                   const std::vector<size_t>& frontier ): _count(count), _frontier(frontier) {}

    void updateWithDeriv( const CountedDerivation& d ) {
        updateCount( d.count );
        updateFrontier( d.bp );
    }

    void updateFrontier( const std::vector<size_t>& bp ) {
        assert(bp.size() == _frontier.size());
        for( size_t ti = 0; ti < _frontier.size(); ++ti ) {
            _frontier[ti] = std::max(bp[ti], _frontier[ti]);
        }
    }

    void updateCount( const cl_I& nc ) { _count += nc; }
    const cl_I& count() { return _count; }
    const std::vector<size_t>& frontier() { return _frontier; }
    void freeDerivations() { _frontier.clear(); }

private:
    std::vector<size_t> _frontier;
    cl_I _count;
};

// lazy

class EdgeDerivInfoLazy{
    // typedefs
public:
    typedef std::vector< std::vector<size_t> > FrontierT;
public:
    EdgeDerivInfoLazy( ) : _count(cl_I(0)), _frontier( vector<vector<size_t>>(0,vector<size_t>()) ) {}
    EdgeDerivInfoLazy( const cl_I& count,
                   const std::vector<size_t>& ibp ): _count(count), _frontier({ibp}) {}

    void updateWithDeriv( const CountedDerivation& d ) {
        updateCount( d.count );
        updateFrontier( d.bp );
    }

    void updateFrontier( const std::vector<size_t>& bp ) {
        _frontier.push_back(bp);
    }

    void updateCount( const cl_I& nc ) { _count += nc; }

    void freeDerivations() { _frontier.clear(); }

    const cl_I& count() { return _count; }
    const std::vector<std::vector<size_t>>& frontier() { return _frontier; }

private:
    std::vector<std::vector<size_t>> _frontier;
    cl_I _count;
};

template <typename EdgeDerivInfoT>
class CostClass {

    typedef size_t edgeID_t;
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

            // lazy
            _usedEdges[e].updateWithDeriv( deriv );

            // eager
            //_usedEdges[e].updateFrontier( deriv.bp );
            //_usedEdges[e].updateCount( numSln );
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
        cl_I edgeCount(0);
        if( hasEdge(eid) ) {
          edgeCount = _usedEdges[eid].count();
        }
        return double_approx(edgeCount / _total);
    }

   /** Return the total number of derivations of this cost class
   */
   cl_I total() { return _total; }

   double cost() { return _cost; }

   /** Return the total number of derivations of this cost class
    */
   cl_I total() const { return _total; }

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
    cl_I _total;
};

#endif // COSTCLASS_HPP
