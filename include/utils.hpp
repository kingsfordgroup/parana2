#ifndef UTILS_HPP
#define UTILS_HPP

#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include <queue>
#include <cstdlib>
#include <vector>

namespace Utils {
    using std::unordered_map;
    using std::unordered_set;
    using std::string;
    using std::tuple;
    using std::cout;
    using std::cerr;
    using std::vector;
    using std::priority_queue;
    using std::tie;

    bool advanceElems( vector<size_t>& ptrs, const vector<size_t>& sizes ) {

        assert( ptrs.size() == sizes.size() );
        // Find the last element
        size_t ctr = sizes.size()-1;
        for ( auto rp = ptrs.rbegin(); rp != ptrs.rend(); ++rp ) {
            if ( *rp < sizes[ctr]-1 ) {
                ++(*rp);
                return true;
            } else {
                if ( rp+1 == ptrs.rend() ) {
                    return false;
                } else {
                    (*rp) = 0;
                }
            }
            --ctr;
        }

    }

    template <typename pqT, typename pqCompT>
    bool appendNext( const vector<size_t>& sizes,
                     vector<pqT>& pq,
                     pqCompT& pqComp,
                     std::function< double(const vector<size_t>&) >& computeScore ) {

        if (pq.empty()) { return false; }

        double score; vector<size_t> inds;
        tie(score,inds) = pq.front();
        std::pop_heap( pq.begin(), pq.end(), pqComp ); pq.pop_back();

        auto inQueue = [&] ( const pqT& e ) -> bool {
            auto eInds = get<1>(e);
            auto vsize = eInds.size();

            for ( auto it = 0; it != pq.size(); ++it ) {
                auto oInds = get<1>( pq[it] );
                if ( eInds == oInds ) {
                    return true;
                }
            }

            return false;
        };

        for (size_t i = 0; i < inds.size(); ++i) {
            vector<size_t> newInds(inds);
            newInds[i]++;
            if ( newInds[i] < sizes[i] ) {
                auto newElem = make_tuple( computeScore(newInds), newInds );
                if (!inQueue(newElem) ) {
                    pq.push_back( newElem );
                    std::push_heap( pq.begin(), pq.end(), pqComp );
                }
            }
        }
        return true;
    }

    class TreeInfo {
    public:
        typedef unordered_map<int, unordered_set<string>> NodeMapT;
        typedef unordered_map<int, unordered_set<int>> NodeIndMapT;
        typedef unordered_map<int, tuple<double,double>> NodeIntervalMapT;
        NodeIndMapT subnodes;
        NodeIndMapT leaves;
        NodeMapT enets;
        NodeIntervalMapT extantInterval;
        string name;

        TreeInfo( const string n ) : name(n), subnodes(NodeIndMapT()),
                                     leaves(NodeIndMapT()),
                                     enets(NodeMapT()),
                                     extantInterval(NodeIntervalMapT()) { cout << "creating tree info: " << name << "\n"; }

        bool inSubnodesOf( int u,  int v ) {
            return subnodes[u].find(v) != subnodes[u].end();
        }
        bool inSubnodesOf( int u,  int v ) const {
            return subnodes.find(u) != subnodes.cend() && subnodes.find(u)->second.find(v) != subnodes.find(u)->second.cend();
        }

        double intervalDistance( int u, int v ) {
            //return static_cast<const TreeInfo&>(*this).intervalDistance(u,v);
            auto birthU = get<0>( extantInterval[u] );
            auto deathU = get<1>( extantInterval[u] );
            auto birthV = get<0>( extantInterval[v] );
            auto deathV = get<1>( extantInterval[v] );
            // If the ranges overlap, the distance is 0
            if ( (birthU <= deathV) && (deathU >= birthV) ) { return 0.0; }
            // otherwise
            if ( birthU >= deathV ) {
                return birthU - deathV;
            }

            if ( birthV >= deathU ) {
                return birthV - deathU;
            }
        }
        /*
        double intervalDistance( int u, int v ) {
            return static_cast<const TreeInfo&>(*this).intervalDistance(u,v);
            }*/
    };


    bool differentExtantNetworks( const TreeInfo& tinfo, int u, int v ) {
        if ( u == v ) { return false; }
        for (auto uit = tinfo.enets.find(u)->second.cbegin(); uit != tinfo.enets.find(u)->second.cend(); ++uit) {
            if ( tinfo.enets.find(v)->second.find(*uit) != tinfo.enets.find(v)->second.cend() ) {
                return false;
            }
        }
        return true;
    }

}

#endif // UTILS_HPP
