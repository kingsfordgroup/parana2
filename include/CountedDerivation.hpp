#ifndef COUNTED_DERIVATION_HPP
#define COUNTED_DERIVATION_HPP

#include <tuple>
#include <vector>
#include <iostream>
#include <unordered_set>
#include <boost/lexical_cast.hpp>
#include <boost/functional/hash.hpp>
#include <boost/multiprecision/gmp.hpp>

#include <cln/rational.h>
#include <cln/integer.h>
#include <cln/float.h>
#include <cln/real.h>
#include <google/dense_hash_set>

#include "ParanaCommon.hpp"

using cln::cl_I;
using boost::multiprecision::mpz_int;
using std::tuple;
using std::unordered_set;
using std::get;
using std::vector;
using std::string;
using std::make_tuple;
using boost::lexical_cast;

class CountedDerivation {
    //typedef mpz_int BigInt;
    typedef cl_I BigInt;
    
public:
    double cost;
    size_t edge;
    vector<size_t, CustomAllocator<size_t>> bp;
    BigInt count;


    CountedDerivation () : cost(0.0), edge(0), bp(vector<size_t, CustomAllocator<size_t>>()), count(BigInt(0)) { }

    CountedDerivation( double _cost, size_t _edge, const vector<size_t, CustomAllocator<size_t>>& _bp, const BigInt& dcount) :
        cost(_cost), edge(_edge), bp(std::move(_bp)), count(dcount) {}

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
    class hash<CountedDerivation> {
    public:
        std::size_t operator()(const CountedDerivation& d) const {
            return d.hashCode();
        }
    };
}

#endif // COUNTED_DERIVATION_HPP
