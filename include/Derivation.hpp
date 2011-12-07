#ifndef DERIVATION_HPP
#define DERIVATION_HPP

#include <tuple>
#include <vector>
#include <unordered_set>
#include <boost/lexical_cast.hpp>
#include <boost/functional/hash.hpp>
#include <google/dense_hash_set>

using google::dense_hash_set;      // namespace where class lives by

using std::tuple;
using std::unordered_set;
using std::get;
using std::vector;

using boost::lexical_cast;

namespace std {
    template<>
    class hash<tuple<int,int,string>> {
    public:
        std::size_t operator()(const tuple<int,int,string>& d) const {
            size_t seed = 0;
            boost::hash_combine(seed, get<0>(d));
            boost::hash_combine(seed, get<1>(d));
            boost::hash_combine(seed, get<2>(d));
            return seed;
        }
    };
}


class Derivation {

public:
    typedef tuple<int,int,string> flipT;
    double cost;
    size_t target;
    vector<size_t> bp;
    dense_hash_set<flipT,std::hash<flipT>> flips;

    friend ostream& operator<<(ostream& output, const Derivation& d) {
        string bpStr;
        for ( auto t : d.bp ) {
            bpStr += lexical_cast<string>(t) + " ";
        }
        string flipStr;
        for ( auto f : d.flips ) {
            flipStr += "[" +lexical_cast<string>(get<0>(f)) + ", " + lexical_cast<string>(get<1>(f)) + " : " + get<2>(f) + "] ";
        }
        output << " cost : " <<  d.cost << ", using edge " << d.target << ", with backpointers " << bpStr << "\n";
        output << " Flips = " << flipStr;

        return output;
    }

    Derivation () : cost(0.0), target(0), bp(vector<size_t>()), flips(dense_hash_set<flipT,std::hash<flipT>>()) {
        flips.set_empty_key(make_tuple(-1,-1,""));
    }

    Derivation( double _cost, size_t _target, const vector<size_t>& _bp, const dense_hash_set<flipT,std::hash<flipT>>& _flips ) :
        cost(_cost), target(_target), bp(_bp), flips(_flips) {}

    size_t hashCode() const {
        size_t seed = 0;
        boost::hash_combine(seed, cost);
        boost::hash_combine(seed, target);
        for ( auto p : bp ) { boost::hash_combine(seed, p); }
        for ( auto f : flips ) {
            boost::hash_combine(seed, get<0>(f));
            boost::hash_combine(seed, get<1>(f));
            boost::hash_combine(seed, get<2>(f));
        }
        return seed;
    }
};

namespace std {
    template<>
    class hash<Derivation> {
    public:
        std::size_t operator()(const Derivation& d) const {
            return d.hashCode();
        }
    };
}


#endif // DERIVATION_HPP
