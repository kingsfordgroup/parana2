#ifndef FWHGRAPH_HPP
#define FWHGRAPH_HPP

#include <vector>
#include <unordered_set>
#include <boost/bimap.hpp>
#include <boost/bimap/unordered_set_of.hpp>
#include <boost/bimap/vector_of.hpp>
#include <functional>
#include <set>
#include <boost/lexical_cast.hpp>
#include "FlipKey.hpp"

using std::vector;
using std::unordered_set;
using boost::bimap;
using boost::bimaps::vector_of;
using boost::bimaps::tagged;
using boost::bimaps::unordered_set_of;
using boost::lexical_cast;
using std::hash;
using std::set;

class HyperEdge {
public:

    friend ostream& operator<<(ostream& output, const HyperEdge& e) {
        string tailStr;
        for ( auto t : e.tail() ) {
            tailStr += lexical_cast<string>(t) + " ";
        }

        output << "[ Head : " << e.head() << " Tail : { " << tailStr << "} ";
        return output;
    }


    HyperEdge( const vector<size_t>& tail, const size_t& head, double weight ) :
        _head(head), _tail(tail), _weight(weight) {}

    const size_t& head() const { return _head; }
    const vector<size_t>& tail() const { return _tail; }
    double weight() const { return _weight; }

    inline bool operator==( const HyperEdge& e ) const {
        if (_head != e._head) { return false;}
        vector<size_t> res; res.reserve(_tail.size());
        std::set_symmetric_difference( e._tail.begin(), e._tail.end(),
                                       _tail.begin(), _tail.end(),
                                       std::back_inserter(res) );
        return res.size() == 0;
    }



private:
    size_t _head;
    vector<size_t> _tail;
    double _weight;
};

/** Make HyperEdges hashable **/
namespace std {
    template<>
    class hash<HyperEdge> {
    public:
        std::size_t operator()(const HyperEdge& e) const {
            hash<size_t> h;
            size_t r = h(e.head());
            size_t m = 1;
            for ( auto t : e.tail() ) {
                r ^= m*h(t); m += 1;
            }
            return r;
        }
    };
}

/** Make HyperEdges hashable **/
namespace boost {
    template<>
    class hash<HyperEdge> {
    public:
        std::size_t operator()(const HyperEdge& e) const {
            hash<size_t> h;
            size_t r = h(e.head());
            size_t m = 1;
            for ( auto t : e.tail() ) {
                r ^= m*h(t); m += 1;
            }
            return r;
        }
    };
}





class ForwardHypergraph {
    typedef bimap<
        vector_of< size_t >,
        unordered_set_of< FlipKey >
        > VertexBiMapT;
    struct edgeT {};
    struct idT {};
    typedef bimap<
        vector_of< tagged<size_t, idT> >,
        unordered_set_of< tagged<HyperEdge, edgeT> >
        > EdgeBiMapT;

    typedef VertexBiMapT::value_type VertexBiKey;
    typedef EdgeBiMapT::value_type EdgeBiMapKey;
    typedef vector< unordered_set<size_t>* > AdjSetT;

public:
    ForwardHypergraph() :
        _vBiMap( VertexBiMapT() ), _eBiMap( EdgeBiMapT() ), _vertices( AdjSetT() ) {}

    bool addVertex( const FlipKey& k ) {
        auto it = _vBiMap.right.find(k);
        if ( it == _vBiMap.right.end() ) {
            //cout << "ADDING " << k << " to the hash\n";
            size_t vid = _vertices.size();
            //cout << "SIZE OF _vBiMap is " << vid << "\n";
            _vBiMap.push_back( VertexBiMapT::value_type(vid, k) );
            _vertices.push_back( new unordered_set<size_t>() );
            return true;
        } // else {//cout << "ALREADY FOUND " << k << " in the hash!\n";}
        return false;
    }

    bool addEdge( const std::vector<FlipKey>& tail, const FlipKey& head, double weight ) {
        auto hind = _vBiMap.right.find(head)->second ;

        vector<size_t> tailInds; tailInds.reserve(tail.size());
        for( auto t : tail ) { tailInds.push_back( _vBiMap.right.find(t)->second ); }
        sort( tailInds.begin(), tailInds.end() );

        HyperEdge edge(tailInds, hind, weight);
        auto it = _eBiMap.right.find(edge);
        if ( it == _eBiMap.right.end() ) {
            size_t eid = _eBiMap.size();
            _eBiMap.push_back( EdgeBiMapT::value_type(eid, edge) );
            (_vertices[hind])->insert( eid );
            //cout << "Adding hyperedge " << edge << "\n";
            return true;
        }
        return false;
    }

    const size_t& getHead( const size_t& eInd ) { return (_eBiMap.left.begin()+eInd)->second.head(); } //.right.head(); }
    const vector<size_t>& getTail( const size_t& eInd) { return (_eBiMap.left.begin()+eInd)->second.tail();}//.right.tail(); }

    const FlipKey& vertex( const size_t& vid ) const {return (_vBiMap.left.begin()+vid)->second;}
    const HyperEdge& edge( const size_t& eid ) const {return (_eBiMap.left.begin()+eid)->second;}

    const size_t& index( const FlipKey& k ) const { return (_vBiMap.right.find(k)->second); }

    size_t order() const { return _vBiMap.size(); }
    size_t size() const { return _eBiMap.size(); }
    const unordered_set<size_t>& incident( const size_t& hind ) {
        if ( hind < _vertices.size() ) {
            return *(_vertices[hind]);
        } else {
            cerr << "ERROR\n";
            exit(1);
        }
    }

private:
    VertexBiMapT _vBiMap;
    EdgeBiMapT _eBiMap;
    AdjSetT _vertices;
};

#endif // FWHGRAPH_HPP
/* -*- c-file-style: "linux" -*- */
