#ifndef COUNTED_DERIVATION_HPP
#define COUNTED_DERIVATION_HPP

#include <tuple>
#include <vector>
#include <iostream>
#include <unordered_set>
#include <boost/lexical_cast.hpp>
#include <boost/functional/hash.hpp>
#include <google/dense_hash_set>

using std::tuple;
using std::unordered_set;
using std::get;
using std::vector;
using std::string;
using std::make_tuple;
using boost::lexical_cast;


class LogCountedDerivation {

public:
    double cost;
    size_t edge;
    vector<size_t> bp;
    double logCount;


    CountedDerivation () : cost(0.0), edge(0), bp(vector<size_t>()), logCount(0.0) { }

    CountedDerivation( double _cost, size_t _edge, const vector<size_t>& _bp, const double& dcount) :
        cost(_cost), edge(_edge), bp(_bp), logCount(dcount) {}

    size_t hashCode() const {
        size_t seed = 0;
        boost::hash_combine(seed, cost);
        boost::hash_combine(seed, edge);
        for ( auto p : bp ) { boost::hash_combine(seed, p); }
        //	boost::hash_combine(seed, count);
        return seed;
    }
};

namespace std {
    template<>
    class hash<LogCountedDerivation> {
    public:
        std::size_t operator()(const LogCountedDerivation& d) const {
            return d.hashCode();
        }
    };
}

#endif // COUNTED_DERIVATION_HPP
