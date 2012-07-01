#ifndef UTILS_HPP
#define UTILS_HPP

#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include <queue>
#include <cstdlib>
#include <vector>

// The logger library
#include "cpplog.hpp"

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

    template <typename pqT, typename pqCompT<pqT>>
    bool appendNext( const vector<size_t>& sizes,
                     vector<pqT>& pq,
                     pqCompT<pqT>& pqComp,
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
                pq.push_back( make_tuple( computeScore(newInds), newInds ) );
                std::push_heap( pq.begin(), pq.end(), pqComp );
            }
            if (inds[i] != 1) {
                break;
            }
         }
            /*
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
            */
        return true;
    }

    template <typename pqT, typename pqCompT<pqT>>
    bool appendNextWithEdge( const vector<size_t>& sizes,
                             vector<pqT>& pq,
                             pqCompT<pqT>& pqComp,
                             std::function< double(const size_t& eid, const vector<size_t>&) >& computeScore ) {

        if (pq.empty()) { return false; }

        double score; size_t eid; vector<size_t> inds;
        tie(score,inds) = pq.front();
        std::pop_heap( pq.begin(), pq.end(), pqComp ); pq.pop_back();

         for (size_t i = 0; i < inds.size(); ++i) {
            vector<size_t> newInds(inds);
            newInds[i]++;
            if ( newInds[i] < sizes[i] ) {
                pq.push_back( make_tuple( computeScore(eid, newInds), eid, newInds ) );
                std::push_heap( pq.begin(), pq.end(), pqComp );
            }
            if (inds[i] != 1) {
                break;
            }
         }
        return true;
    }


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
