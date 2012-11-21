#include "MultiOpt.hpp"
#include "TreeUtils.hpp"
#include "ProgressDisplay.hpp"
#include "model.hpp"

#include <stack>
#include <cmath>
#include <utility>
#include <boost/timer/timer.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <cln/float.h>

/** Google's dense hash set and hash map **/
#include <google/dense_hash_set>
#include <google/dense_hash_map>
#include <google/sparse_hash_map>

// For logging
#include "cpplog.hpp"

namespace MultiOpt {
using Utils::round3;

const auto none = make_tuple(false, false);
const auto fw = make_tuple(true, false);
const auto rev = make_tuple(false, true);
const auto both = make_tuple(true, true);

static unordered_map<tuple<bool, bool>, string> flipStrMap = {
    {none, "n"},
    {fw, "f"},
    {rev, "r"},
    {both, "b"}
};

static flipMapT flipDict = {
    {none, { {both, "b+"}, {fw, "f+"}, {rev, "r+"}, {none, "n"} } },
    {fw,   { {both, "r+"}, {rev, "f-r+"}, {none, "f-"}, {fw, "n"} } },
    {rev,  { {both, "f+"}, {rev, "n"}, {none, "r-"}, {fw, "f+r-"} } },
    {both, { {both, "n"}, {fw, "r-"}, {rev, "f-"}, {none, "b-"} } }
};

costMapT getCostDict ( double cc, double dc, bool directed ) {
    auto undirected = ! directed;
    auto none = make_tuple(false, false);
    auto fw = make_tuple(true, false);
    auto rev = make_tuple(false, true);
    auto both = make_tuple(true, true);

    costMapT costDict;
    costDict[ none ][ none ] = make_tuple( 0.0, "n" );
    costDict[ none ][ fw ] = make_tuple( cc, "f+" );
    costDict[ none ][ rev ] = make_tuple( cc, "r+" );
    costDict[ none ][ both ] = make_tuple( undirected ? cc : 2 * cc , "b+" );

    costDict[ rev ][ none ] = make_tuple( dc, "r-" );
    costDict[ rev ][ fw ] = make_tuple( cc + dc, "f+r-" );
    costDict[ rev ][ rev ] = make_tuple( 0.0, "n" );
    costDict[ rev ][ both ] = make_tuple( cc , "f+" );

    costDict[ fw ][ none ] = make_tuple( dc, "f-" );
    costDict[ fw ][ fw ] = make_tuple( 0.0, "n" );
    costDict[ fw ][ rev ] = make_tuple( cc + dc, "f-r+" );
    costDict[ fw ][ both ] = make_tuple( cc , "r+" );

    costDict[ both ][ none ] = make_tuple( undirected ? dc : 2 * dc, "b-" );
    costDict[ both ][ fw ] = make_tuple( dc, "r-" );
    costDict[ both ][ rev ] = make_tuple( dc, "f-" );
    costDict[ both ][ both ] = make_tuple( 0.0, "n");

    return costDict;
}

costMapFunT getCostFunDict ( double cc, double dc, bool directed ) {
    auto undirected = ! directed;
    auto none = make_tuple(false, false);
    auto fw = make_tuple(true, false);
    auto rev = make_tuple(false, true);
    auto both = make_tuple(true, true);

    costMapFunT costDict;
    costDict[ none ][ none ] = [ = ](const double & pf, const double & pr) {
        return make_tuple( 0.0, "n" );
    };
    costDict[ none ][ fw ] = [ = ](const double & pf, const double & pr) {
        return make_tuple( pf * cc, "f+" );
    };
    costDict[ none ][ rev ] = [ = ](const double & pf, const double & pr) {
        return make_tuple( pr * cc, "r+" );
    };
    costDict[ none ][ both ] = [ = ](const double & pf, const double & pr) {
        return make_tuple( undirected ? pf * cc : (pf + pr) * cc, "b+" );
    };


    costDict[ rev ][ none ] = [ = ](const double & pf, const double & pr) {
        return make_tuple( pr * dc, "r-" );
    };
    costDict[ rev ][ fw ] = [ = ](const double & pf, const double & pr) {
        return make_tuple( pf * cc + pr * dc, "f+r-" );
    };
    costDict[ rev ][ rev ] = [ = ](const double & pf, const double & pr) {
        return make_tuple( (1.0 - pr) * dc, "n");
    };
    costDict[ rev ][ both ] = [ = ](const double & pf, const double & pr) {
        return make_tuple( pf * cc + (1.0 - pr) * dc , "f+");
    };

    costDict[ fw ][ none ] = [ = ](const double & pf, const double & pr) {
        return make_tuple( pf * dc, "f-");
    };
    costDict[ fw ][ fw ] = [ = ](const double & pf, const double & pr) {
        return make_tuple( (1.0 - pf) * dc, "n");
    };
    costDict[ fw ][ rev ] = [ = ](const double & pf, const double & pr) {
        return make_tuple( pf * dc + pr * cc, "f-r+");
    };
    costDict[ fw ][ both ] = [ = ](const double & pf, const double & pr) {
        return make_tuple( (1.0 - pf) * dc + pr * cc, "r+");
    };

    costDict[ both ][ none ] = [ = ](const double & pf, const double & pr) {
        return make_tuple( undirected ? dc : (pf + pr) * dc, "b-" );
    };
    costDict[ both ][ fw ] = [ = ](const double & pf, const double & pr) {
        return make_tuple( (1.0 - pf) * dc + pr * dc, "r-" );
    };
    costDict[ both ][ rev ] = [ = ](const double & pf, const double & pr) {
        return make_tuple( pf * dc + (1.0 - pr) * dc, "f-" );
    };
    costDict[ both ][ both ] = [ = ](const double & pf, const double & pr) {
        return make_tuple( undirected ? (1.0 - pf) * dc : (1.0 - pf) * dc + (1.0 - pr) * dc, "n" );
    };

    return costDict;
}

selfCostMapT getSelfLoopCostDict ( double cc, double dc, bool directed ) {
    auto undirected = ! directed;
    auto none = make_tuple(false, false);
    auto fw = make_tuple(true, false);
    auto rev = make_tuple(false, true);
    auto both = make_tuple(true, true);

    selfCostMapT selfLoopCostDict;
    selfLoopCostDict[ none ][ true ] = make_tuple( cc, "b+" );
    selfLoopCostDict[ none ][ false ] = make_tuple( 0.0, "n" );

    selfLoopCostDict[ rev ][ true ] = make_tuple( 0.0, "n" );
    selfLoopCostDict[ rev ][ false ] = make_tuple( dc, "b-" );

    selfLoopCostDict[ fw ][ true ] = make_tuple( 0.0, "n" );
    selfLoopCostDict[ fw ][ false ] = make_tuple( dc, "b-" );

    selfLoopCostDict[ both ][ true ] = make_tuple( 0.0, "n" );
    selfLoopCostDict[ both ][ false ] = make_tuple( dc, "b-" );
    return selfLoopCostDict;
}


selfCostMapFunT getSelfLoopCostFunDict ( double cc, double dc, bool directed ) {
    auto undirected = ! directed;
    auto none = make_tuple(false, false);
    auto fw = make_tuple(true, false);
    auto rev = make_tuple(false, true);
    auto both = make_tuple(true, true);

    selfCostMapFunT selfLoopCostDict;
    selfLoopCostDict[ none ][ true ] = [ = ] (const double & p) {
        return make_tuple( p * cc, "b+" );
    };
    selfLoopCostDict[ none ][ false ] = [ = ] (const double & p) {
        return make_tuple( 0.0, "n" );
    };

    selfLoopCostDict[ rev ][ true ] = [ = ] (const double & p) {
        return make_tuple( (1.0 - p) * dc, "n" );
    };
    selfLoopCostDict[ rev ][ false ] = [ = ] (const double & p) {
        return make_tuple( p * dc, "b-" );
    };

    selfLoopCostDict[ fw ][ true ] = [ = ] (const double & p) {
        return make_tuple( (1.0 - p) * dc, "n" );
    };
    selfLoopCostDict[ fw ][ false ] = [ = ] (const double & p) {
        return make_tuple( p * dc, "b-" );
    };

    selfLoopCostDict[ both ][ true ] = [ = ] (const double & p) {
        return make_tuple( (1.0 - p) * dc, "n" );
    };
    selfLoopCostDict[ both ][ false ] = [ = ] (const double & p) {
        return make_tuple( dc, "b-" );
    }; // was p*dc
    return selfLoopCostDict;
}

string flipType( const FlipKey &hvert, const FlipKey &tvert ) {
    return flipDict[ make_tuple(hvert.f(), hvert.r())][ make_tuple(tvert.f(), tvert.r()) ];
}

template<typename GT>
std::unordered_set<int> projectToReversedGraph( unique_ptr<ForwardHypergraph> &H, GT &G ) {
    auto M = H->size();
    std::unordered_set<int> vset;
    //for ( size_t vi = 0; vi < M; ++vi ) { vset.insert(vi); }

    for ( size_t eid = 0; eid < M; ++eid ) {
        auto e = H->edge(eid);
        auto h = e.head();
        vset.insert(h);
        for ( auto t : e.tail() ) {
            vset.insert(t);
            add_edge(h, t, G);
        }
    }
    //cerr << "G size = " << num_edges(G) << ", G order = " << num_vertices(G) << "\n";
    //cerr << "Inserted " << vset.size() << " vertices \n";
    return vset;
    /*cerr << "DIFF = {";
    for (const auto& c : vset) { cerr << c << ", "; }
    cerr << "}\n";
    */
}

void topologicalOrder( unique_ptr<ForwardHypergraph> &H, vector<size_t> &order ) {
    using boost::adjacency_list;
    using boost::vecS;
    using boost::directedS;

    cpplog::FileLogger log( "log.txt", true );

    typedef adjacency_list<vecS, vecS, directedS> graphT;
    graphT G( H->order() );

    vector<size_t> torder;
    std::unordered_set<int> vset = projectToReversedGraph( H, G );
    boost::topological_sort(G, std::back_inserter(torder));

    // Remove all of the vertices from "order" if they are not in vset
    LOG_INFO(log) << "ORDER BEFORE = " << torder.size() << "\n";
    std::copy_if(torder.begin(), torder.end(), std::back_inserter(order), [ = ](const int & id) -> bool { return vset.find(id) != vset.end(); } );
    //std::copy()
    LOG_INFO(log) << "ORDER AFTER = " << order.size() << "\n";
}

void topologicalOrderHelper( unique_ptr<ForwardHypergraph> &H, size_t rootInd, 
                             std::unordered_set<size_t> seen, vector<size_t> &order ) {
    if ( seen.find(rootInd) == seen.end() ) {
        if( H->isLeaf(rootInd) ) {
            order.push_back(rootInd);
            seen.insert(rootInd);
        } else {
            // Visit all of the children
            for ( auto e : H->incident(rootInd) ) {
                auto edge = H->edge(e);
                for ( auto v : edge.tail() ) {
                    topologicalOrderHelper(H, v, seen, order);
                }
            }
            // Then visit and mark the parent
            order.push_back(rootInd);
            seen.insert(rootInd);
        }
    }
}

void topologicalOrderQueue( unique_ptr<ForwardHypergraph> &H, size_t rootInd, vector<size_t> &order ) {
    std::unordered_set<size_t> vset;
    topologicalOrderHelper(H, rootInd, vset, order);
    if ( vset.size() < H->order() ) {
        for( size_t i = 0; i < H->order(); ++i ) {
            if ( vset.find(i) == vset.end() ) {
                std::cerr << "missed " << H->vertex(i) << "\n";
            }
        }
    }
}

/**
 *  Compute the penalty for this edge to exist based on difference
 *  between the existence intervals of the endpoints and the
 *  penalty factor.
 */
template <typename T>
double existencePenalty( const TreeInfo &ti, const T &vert, const double &penalty, const double &travWeight ) {

    auto dist = ti.intervalDistance( vert.u(), vert.v() );
    return (travWeight > 0.0) ?  dist / std::exp(-dist * penalty) : 0.0;
    //return (travWeight > 0.0 && dist > 0.0 ) ?  (dist*penalty) : 0.0; //1.0 - std::exp(-penalty*dist) : 0.0;
    //return (travWeight > 0.0) ? ( (penalty*dist) / (1.0 + std::exp(-dist*penalty))) : 0.0;
    //return (travWeight > 0.0) ? 1.0 + ( dist*penalty) : 1.0;
    //return (dist > 0.0 && travWeight > 0.0 ) ? (3.0 + (dist*penalty)) : 1.0;// + (penalty * dist);

}


unique_ptr<ForwardHypergraph>  buildMLSolutionSpaceGraph( const TreePtrT &t,
        const TreeInfo &ti,
        Model &model,
        bool directed ) {

    cpplog::FileLogger log( "log.txt", true );

    boost::timer::auto_cpu_timer timer;

    Google<int>::Set leafSet;
    leafSet.set_empty_key(-1);
    for ( auto l : t->getLeavesId() ) {
        leafSet.insert(l);
    }
    auto isLeaf = [&]( const int & nid ) -> bool { return leafSet.find(nid) != leafSet.end(); }; //t->isLeaf(nid); };
    auto isInternal = [&]( const int & nid ) -> bool { return leafSet.find(nid) == leafSet.end(); };
    auto isLost = [&]( int nid ) -> bool { return (t->getNodeName(nid)).find("LOST") != std::string::npos; };

    auto rootId = t->getRootId();
    auto rootName = t->getNodeName( rootId );

    auto fauxRoot = "#preroot#";

    unique_ptr<ForwardHypergraph> slnSpaceGraph( new ForwardHypergraph() );

    //costMapT costMap( getCostDict(cc, dc, directed) );
    //selfCostMapT selfLoopCostMap( getSelfLoopCostDict(cc, dc, directed) );

    auto addIncomingHyperedge = [ = , &slnSpaceGraph, &t] (const FlipKey & k ) {
        typedef std::vector<FlipKey> FlipKeyVec;

        auto addCartesianProduct = [&] ( std::vector<FlipKeyVec> &iterKeys, std::vector<FlipKey> &fixedKeys ) -> void {
            std::vector<size_t> indexes( iterKeys.size(), 0 );
            bool done = false;//iterKeys.size() == 0;
            while ( !done ) {
                FlipKeyVec res;
                for ( size_t i = 0; i < indexes.size(); ++i ) {
                    res.push_back(iterKeys[i][indexes[i]]);
                }
                res.insert( res.end(), fixedKeys.begin(), fixedKeys.end() );
                double w = 1.0;
                if ( res.size() > 0 ) {
                    slnSpaceGraph->addEdge( res, k, w );
                }
                if ( indexes.size() == 0 ) {
                    return;
                }
                int j = indexes.size() - 1;
                while ( true ) {
                    indexes[j] += 1;
                    if ( indexes[j] < iterKeys[j].size() ) {
                        break;
                    }
                    indexes[j] = 0;
                    j -= 1;
                    if ( j < 0 ) {
                        return;
                    }
                }
            }
        };

        auto dirKey = k.getDirTuple();
        double canonicalDerivCost = 0.0;
        const int LN = k.u();
        const int RN = k.v();

        // Self loop
        if (k.arity() == 1) {
            auto selfNode = LN;
            if (isInternal(selfNode)) {
               // Get the children nodes
                int LLN = t->getSonsId(selfNode)[0];
                int LRN = t->getSonsId(selfNode)[1];
                if ( LLN > LRN ) {
                    std::swap(LLN, LRN);
                }

                // Possible states
                // {LLN,LLN}  {LRN, LRN}  {LLN, LRN}
                //     0          0           0
                //     0          0           1
                //     0          1           0
                //     0          1           1
                //     1          0           0
                //     1          0           1
                //     1          1           0
                //     1          1           1
                std::vector<FlipKey> leftKeys;
                std::vector<FlipKey> rightKeys;
                std::vector<FlipKey> betweenKeys;

                if ( !isLost(LLN) ) {
                    leftKeys.push_back( {LLN, LLN, true, true});
                    leftKeys.push_back( {LLN, LLN, false, false});
                }
                if ( !isLost(LRN) ) {
                    rightKeys.push_back( {LRN, LRN, true, true});
                    rightKeys.push_back( {LRN, LRN, false, false});

                }
                if (! differentExtantNetworks(ti, LLN, LRN) && ! (isLost(LLN) || isLost(LRN)) ) {
                    betweenKeys.push_back( {LLN, LRN, true, true});
                    betweenKeys.push_back( {LLN, LRN, false, false});
                }

                std::vector< FlipKeyVec > iterKeys;
                if ( leftKeys.size() > 0 ) {
                    iterKeys.push_back(leftKeys);
                }
                if ( rightKeys.size() > 0 ) {
                    iterKeys.push_back(rightKeys);
                }
                if ( betweenKeys.size() > 0 ) {
                    iterKeys.push_back( betweenKeys);
                }

                std::vector<FlipKey> emptyVector;
                addCartesianProduct(iterKeys, emptyVector);
            }
        } else {
            std::vector<FlipKey> llnChoices;
            std::vector<FlipKey> lrnChoices;
            std::vector<FlipKey> rlnChoices;
            std::vector<FlipKey> rrnChoices;
            if ( isInternal(LN) ) {
                // Childrent of left node
                int LLN = t->getSonsId(LN)[0];
                int LRN = t->getSonsId(LN)[1];
                if ( LLN > LRN ) {
                    std::swap(LLN, LRN);
                }

                if ( !differentExtantNetworks(ti, LLN, RN) && !(isLost(LLN) || isLost(RN)) ) {
                    llnChoices.push_back( FlipKey(LLN, RN, true, true) );
                    llnChoices.push_back( FlipKey(LLN, RN, false, false) );
                }

                if ( !differentExtantNetworks(ti, LRN, RN) && !(isLost(LRN) || isLost(RN)) ) {
                    lrnChoices.push_back( FlipKey(LRN, RN, true, true) );
                    lrnChoices.push_back( FlipKey(LRN, RN, false, false) );
                }

            }

            if ( isInternal(RN) ) {
                // Childrent of right node
                int RLN = t->getSonsId(RN)[0];
                int RRN = t->getSonsId(RN)[1];
                if ( RLN > RRN ) {
                    std::swap(RLN, RRN);
                }

                if ( !differentExtantNetworks(ti, RLN, LN) && !(isLost(RLN) || isLost(LN)) ) {
                    rlnChoices.push_back( FlipKey(RLN, LN, true, true) );
                    rlnChoices.push_back( FlipKey(RLN, LN, false, false) );
                }

                if ( !differentExtantNetworks(ti, RRN, LN) && !(isLost(RRN) || isLost(LN)) ) {
                    rrnChoices.push_back( FlipKey(RRN, LN, true, true) );
                    rrnChoices.push_back( FlipKey(RRN, LN, false, false) );
                }

            }


            std::vector< FlipKeyVec > iterKeys;
            std::vector< FlipKey > fixedKeys;

            /**  Covers the following cases
            *    T,T,F,F
            *    T,F,F,F
            *    F,T,F,F
            *    F,F,F,F
            */
            if ( llnChoices.size() > 0 ) {
                iterKeys.push_back(llnChoices);
            }
            if ( lrnChoices.size() > 0 ) {
                iterKeys.push_back(lrnChoices);
            }
            if ( rlnChoices.size() > 0 ) {
                fixedKeys.push_back(rlnChoices[1]);
            }
            if ( rrnChoices.size() > 0 ) {
                fixedKeys.push_back(rrnChoices[1]);
            }
            std::vector< FlipKey > none;
            addCartesianProduct(iterKeys, none);

            iterKeys.clear(); fixedKeys.clear();
            if ( llnChoices.size() > 0 ) {
                fixedKeys.push_back(llnChoices[1]);
            }
            if ( lrnChoices.size() > 0 ) {
                fixedKeys.push_back(lrnChoices[1]);
            }
            /** Covers the following cases
            *   F,F,T,T
            *   F,F,F,T
            *   F,F,T,F
            */
            if ( rlnChoices.size() > 0 and rrnChoices.size() > 0 ) {
                std::vector<FlipKey> res;
                res = {rlnChoices[0], rrnChoices[0]};
                //res.insert(res.end(), fixedKeys.begin(), fixedKeys.end());
                if ( res.size() > 0 ) {
                    slnSpaceGraph->addEdge( res, k, 1.0);
                }

                res = {rlnChoices[1], rrnChoices[0]};
                //res.insert(res.end(), fixedKeys.begin(), fixedKeys.end());
                if ( res.size() > 0 ) {
                    slnSpaceGraph->addEdge( res, k, 1.0);
                }

                res = {rlnChoices[0], rrnChoices[1]};
                //res.insert(res.end(), fixedKeys.begin(), fixedKeys.end());
                if ( res.size() > 0 ) {
                    slnSpaceGraph->addEdge( res, k, 1.0);
                }

                res = {rlnChoices[1], rrnChoices[1]};
                if ( res.size() > 0 ) {
                    slnSpaceGraph->addEdge( res, k , 1.0);
                }
            } else if ( rlnChoices.size() > 0 ) {
                std::vector<FlipKey> res;
                res = {rlnChoices[0]};
                //res.insert(res.end(), fixedKeys.begin(), fixedKeys.end());
                if ( res.size() > 0 ) {
                    slnSpaceGraph->addEdge( res, k, 1.0);
                }
                res = {rlnChoices[1]};
                //res.insert(res.end(), fixedKeys.begin(), fixedKeys.end());
                if ( res.size() > 0 ) {
                    slnSpaceGraph->addEdge( res, k, 1.0);
                }
            } else if ( rrnChoices.size() > 0 ) {
                std::vector<FlipKey> res;
                res = {rrnChoices[0]};
                //res.insert(res.end(), fixedKeys.begin(), fixedKeys.end());
                if ( res.size() > 0 ) {
                    slnSpaceGraph->addEdge( res, k, 1.0);
                }
                res = {rrnChoices[1]};
                //res.insert(res.end(), fixedKeys.begin(), fixedKeys.end());
                if ( res.size() > 0 ) {
                    slnSpaceGraph->addEdge( res, k, 1.0);
                }
            }
        }
    };

    // Add the nodes to the hypergraph

    vector<int> nodes;
    for ( auto n : t->getNodesId() ) {
        nodes.push_back(n);
    }

    auto tbegin = nodes.cbegin();
    auto tend = nodes.cend();
    vector<int>::const_iterator uit = tbegin;
    vector<int>::const_iterator vit = tbegin;

    // For every pair of nodes in the tree
    for ( uit = tbegin; uit != tend; uit++ ) {
        int u = *uit;
        for ( vit = uit; vit != tend; vit++ ) {
            int v = *vit;

            // Add the pair if they meet the following conditions
            if ( (u == v) || // they are the same node (self-loop)
                    (!differentExtantNetworks(ti, u, v) && // they have progeny in the same species
                     !(ti.inSubnodesOf(u, v) ||  ti.inSubnodesOf(v, u)) && // neither is a progeny of the other
                     !(isLost(u) || isLost(v)) )) { // neither of them is lost
                // Add a vertex representing no edge between them
                slnSpaceGraph->addVertex( FlipKey( u, v, false, false ) );
                if ( ! t->isRoot(v) ) {
                    // Add a vertex representing an edge between them (or both edges in a dir. graph)
                    slnSpaceGraph->addVertex( FlipKey( u, v, true, true ) );
                }
                if ( directed && u != v ) {
                    // If this is a directed reconstruction, add both the forward and reverse edges
                    slnSpaceGraph->addVertex( FlipKey( u, v, true, false ) );
                    slnSpaceGraph->addVertex( FlipKey( u, v, false, true ) );
                }

            } // end node addition
        } // end vit
    } // end uit

    // The hypergraph has N vertices
    auto N = slnSpaceGraph->order();
    ProgressDisplay showProgress(N);
    // For each vertex
    for ( size_t i = 0; i < N; ++i, ++showProgress ) {
        auto k = slnSpaceGraph->vertex(i);
        addIncomingHyperedge( k );
    }

    //cerr << "\n";
    LOG_INFO(log) << "Hypergraph size = " << slnSpaceGraph->size() << "\n";
    return slnSpaceGraph;
}


unique_ptr<ForwardHypergraph>  buildSolutionSpaceGraph( const TreePtrT &t,
        const TreeInfo &ti,
        double cc,
        double dc,
        double penalty,
        bool directed ) {

    cpplog::FileLogger log( "log.txt", true );

    boost::timer::auto_cpu_timer timer;

    Google<int>::Set leafSet;
    leafSet.set_empty_key(-1);
    for ( auto l : t->getLeavesId() ) {
        leafSet.insert(l);
    }
    auto isLeaf = [&]( const int & nid ) -> bool { return leafSet.find(nid) != leafSet.end(); }; //t->isLeaf(nid); };
    auto isInternal = [&]( const int & nid ) -> bool { return leafSet.find(nid) == leafSet.end(); };
    auto isLost = [&]( int nid ) -> bool { return (t->getNodeName(nid)).find("LOST") != std::string::npos; };

    auto rootId = t->getRootId();
    auto rootName = t->getNodeName( rootId );

    auto fauxRoot = "#preroot#";

    unique_ptr<ForwardHypergraph> slnSpaceGraph( new ForwardHypergraph() );

    costMapT costMap( getCostDict(cc, dc, directed) );
    selfCostMapT selfLoopCostMap( getSelfLoopCostDict(cc, dc, directed) );

    auto addIncomingHyperedge = [ = , &costMap, &selfLoopCostMap, &slnSpaceGraph, &penalty, &t] (const FlipKey & k, const int & rnode, const int & onode ) {
        auto dirKey = k.getDirTuple();
        double canonicalDerivCost = 0.0;

        // Self loop
        if (k.arity() == 1) {
            auto selfNode = rnode;
            if (isInternal(selfNode)) {
                // Get the children nodes
                int LRN = -1;
                int RRN = -1;
                if (t->getNodeName(t->getSonsId(selfNode)[0]) <  t->getNodeName(t->getSonsId(selfNode)[1])) {
                    LRN = t->getSonsId(selfNode)[0];
                    RRN = t->getSonsId(selfNode)[1];
                } else {
                    LRN = t->getSonsId(selfNode)[1];
                    RRN = t->getSonsId(selfNode)[0];
                }
                // This vertex has 2 incoming hyper edges
                // 1 -- we don't flip the self-loop
                auto noFlipLL = FlipKey( LRN, LRN, k.f(), k.r() );
                auto noFlipRR = FlipKey( RRN, RRN, k.f(), k.r() );
                auto noFlipLR = FlipKey( LRN, RRN, k.f(), k.r() );
                // 2 -- we flip the self loop
                auto dualFlipLL = flipBoth( noFlipLL );
                auto dualFlipRR = flipBoth( noFlipRR );
                auto dualFlipLR = flipBoth( noFlipLR );
                vector<FlipKey> noFlipEdge = vector<FlipKey>();
                vector<FlipKey> dualFlipEdge = vector<FlipKey>();

                if ( !isLost(LRN) ) {
                    noFlipEdge.push_back( noFlipLL );
                    dualFlipEdge.push_back( dualFlipLL );
                }
                if ( !isLost(RRN) ) {
                    noFlipEdge.push_back( noFlipRR );
                    dualFlipEdge.push_back( dualFlipRR );
                }

                //vector<FlipKey> noFlipEdge = { noFlipLL, noFlipRR };
                //vector<FlipKey> dualFlipEdge = { dualFlipLL, dualFlipRR };

                if (! differentExtantNetworks(ti, LRN, RRN) && ! (isLost(LRN) || isLost(RRN)) ) {
                    noFlipEdge.push_back( noFlipLR );
                    dualFlipEdge.push_back( dualFlipLR );
                }
                if ( noFlipEdge.size() > 0 ) {
                    slnSpaceGraph->addEdge( noFlipEdge, k, 0.0 );
                }
                if ( dualFlipEdge.size() > 0 ) {
                    auto w = round3(get<0>(selfLoopCostMap[ make_tuple(k.f(), k.r()) ][ dualFlipLL.f() || dualFlipLL.r() ]));
                    w += existencePenalty(ti, k, penalty, w);
                    if ( std::isfinite(w) ) {
                        slnSpaceGraph->addEdge( dualFlipEdge, k, w );
                    }
                }
            }
        } else {
            if ( isInternal(rnode) ) {
                //cout << t->getNodeName(rnode) << " is INTERNAL\n";
                // Get the children nodes
                //int LRN = t->getSonsId(rnode)[0]; int RRN = t->getSonsId(rnode)[1];
                int LRN = -1;
                int RRN = -1;
                if (t->getNodeName(t->getSonsId(rnode)[0]) <  t->getNodeName(t->getSonsId(rnode)[1])) {
                    LRN = t->getSonsId(rnode)[0];
                    RRN = t->getSonsId(rnode)[1];
                } else {
                    LRN = t->getSonsId(rnode)[1];
                    RRN = t->getSonsId(rnode)[0];
                }
                // This vertex has 2 incoming hyper edges
                // 1 -- we don't flip the self-loop
                auto noFlipL = FlipKey( LRN, onode, k.f(), k.r() );
                auto noFlipR = FlipKey( RRN, onode, k.f(), k.r() );
                // 2 -- we flip the self loop
                auto dualFlipL = flipBoth( noFlipL );
                auto dualFlipR = flipBoth( noFlipR );

                vector<FlipKey> noFlip;
                vector<FlipKey> dualFlip;

                if ( !differentExtantNetworks(ti, LRN, onode) && !(isLost(LRN) || isLost(onode)) ) {
                    noFlip.push_back( noFlipL );
                    dualFlip.push_back( dualFlipL );
                }
                if ( !differentExtantNetworks(ti, RRN, onode) && !(isLost(RRN) || isLost(onode)) ) {
                    noFlip.push_back( noFlipR );
                    dualFlip.push_back( dualFlipR );
                }

                if ( noFlip.size() > 0 ) {
                    slnSpaceGraph->addEdge( noFlip, k, canonicalDerivCost);//0.0 );
                }
                if ( dualFlip.size() > 0 ) {
                    auto w = round3(get<0>(costMap[ make_tuple(k.f(), k.r()) ][ make_tuple(dualFlipL.f(), dualFlipL.r()) ])); //directed ? 2.0 : 1.0;//2.0 ? directed : 1.0;
                    w += existencePenalty(ti, k, penalty, w);
                    if ( std::isfinite(w) ) {
                        slnSpaceGraph->addEdge( dualFlip, k, w );
                    }
                }

                if ( directed ) {
                    auto fwFlipL = flipForward(noFlipL);
                    auto fwFlipR = flipForward(noFlipR);
                    auto revFlipL = flipReverse(noFlipL);
                    auto revFlipR = flipReverse(noFlipR);

                    vector<FlipKey> fwFlip;
                    vector<FlipKey> revFlip;

                    if ( !differentExtantNetworks(ti, LRN, onode) ) {
                        fwFlip.push_back( fwFlipL );
                        revFlip.push_back( revFlipL );
                    }

                    if ( !differentExtantNetworks(ti, RRN, onode) ) {
                        fwFlip.push_back( fwFlipL );
                        revFlip.push_back( revFlipR );
                    }

                    if ( fwFlip.size() > 0 ) {
                        auto w = round3(get<0>(costMap[ make_tuple(k.f(), k.r() ) ][ make_tuple(fwFlipL.f(), fwFlipL.r()) ]));
                        w += existencePenalty(ti, k, penalty, w);
                        if ( std::isfinite(w) ) {
                            slnSpaceGraph->addEdge( fwFlip, k , w );
                        }
                    }
                    if ( revFlip.size() > 0 ) {
                        auto w = round3(get<0>(costMap[ make_tuple(k.f(), k.r() ) ][ make_tuple(revFlipL.f(), revFlipL.r()) ]));
                        w += existencePenalty(ti, k, penalty, w);
                        if ( std::isfinite(w) ) {
                            slnSpaceGraph->addEdge( revFlip, k, w );
                        }
                    }
                }

            }
        }
    };

    // Add the nodes to the hypergraph

    vector<int> nodes;
    for ( auto n : t->getNodesId() ) {
        nodes.push_back(n);
    }

    auto tbegin = nodes.cbegin();
    auto tend = nodes.cend();
    vector<int>::const_iterator uit = tbegin;
    vector<int>::const_iterator vit = tbegin;

    size_t ctr1 = 0;
    size_t ctr2 = 0;
    size_t ctr3 = 0;
    for ( uit = tbegin; uit != tend; uit++ ) {
        int u = *uit;
        for ( vit = uit; vit != tend; vit++ ) {
            int v = *vit;

            if ( (u == v) ||
                    (!differentExtantNetworks(ti, u, v) &&
                     !(ti.inSubnodesOf(u, v) ||  ti.inSubnodesOf(v, u)) &&
                     !(isLost(u) || isLost(v)) )) {
                ctr1 += 1;
                slnSpaceGraph->addVertex( FlipKey( u, v, false, false ) );
                if ( ! t->isRoot(v) ) {
                    ctr2 += 1;
                    slnSpaceGraph->addVertex( FlipKey( u, v, true, true ) );
                }
                if ( directed && u != v ) {
                    ctr3 += 1;
                    slnSpaceGraph->addVertex( FlipKey( u, v, true, false ) );
                    slnSpaceGraph->addVertex( FlipKey( u, v, false, true ) );
                }

            }
        }
    }

    auto N = slnSpaceGraph->order();
    ProgressDisplay showProgress(N);
    //cerr << "Hypergraph size: " << slnSpaceGraph->size() << ", order: " << slnSpaceGraph->order() << "\n";
    for ( size_t i = 0; i < N; ++i, ++showProgress ) {
        /*
        if ( !(i % 1000) || i == N-1 ) {
            cerr << "\r\rProcessed " << i << "/" << N << " nodes";
        }
        */
        auto k = slnSpaceGraph->vertex(i);
        auto u = k.u();
        auto v = k.v();
        auto f = k.f();
        auto r = k.r();
        addIncomingHyperedge( k, u, v );
        if ( k.arity() > 1 ) {
            addIncomingHyperedge( k, v, u );
        }
    }

    /*
    auto fname = "HYPERGRAPH.txt";
    std::fstream out( fname, std::fstream::out | std::fstream::trunc );
    auto nedge = slnSpaceGraph->size();
    // write # of verts
    out << nedge << "\n";
    for( size_t i = 0; i < nedge; ++i ) {
        // for each vertex
        auto edge = slnSpaceGraph->edge(i);
        auto hvert = slnSpaceGraph->vertex(edge.head());
        // write vertex key and # of incident edges
        out << vstr(hvert) << "\t" << edge.weight() << "\t" << edge.tail().size() << "\n";
        for ( auto tind : edge.tail() ) {
            out << vstr( slnSpaceGraph->vertex(tind) ) << "\n";
        }
    }
    out.close();
    */
    //cerr << "\n";
    LOG_INFO(log) << "Hypergraph size = " << slnSpaceGraph->size() << "\n";
    return slnSpaceGraph;
}

template< typename GT >
void MLLeafCostDict( unique_ptr<ForwardHypergraph> &H, TreePtrT &T, GT &G, bool directed, double cc, double dc, slnDictT &slnDict ) {
    /*
      Given the duplication tree T, the root vertex rv, the extant graph G and
      the constraints, fill in the values for the leaf nodes of the hypergraph
    */
    typedef typename boost::graph_traits< GT >::edge_descriptor EdgeT;
    boost::timer::auto_cpu_timer timer;
    auto undirected = !directed;
    auto isLost = [&]( int nid ) -> bool { return (T->getNodeName(nid)).find("LOST") != std::string::npos; };
    // Cost of going from the state encoded by a node to the state of a true graph edge

    auto none = make_tuple(false, false);
    auto fw = make_tuple(true, false);
    auto rev = make_tuple(false, true);
    auto both = make_tuple(true, true);

    auto costDict = getCostDict(cc, dc, directed);
    auto selfLoopCostDict = getSelfLoopCostDict(cc, dc, directed);
    auto costFunDict = getCostFunDict(cc, dc, directed);
    auto selfLoopCostFunDict = getSelfLoopCostFunDict(cc, dc, directed );

    auto N = H->order();
    auto M = H->size();

    // The list of all hypernodes with no descendants
    // We'll have either N^2 or (N^2)/2 leaf hypernodes (depending
    // on directedness)
    auto numExtantNodes = std::distance( vertices(G).first, vertices(G).second );
    auto numConn = ( numExtantNodes * numExtantNodes );
    if (undirected) {
        numConn /= 2;
    }
    vector<size_t> leafHypernodes;
    leafHypernodes.reserve(numConn);

    // For every hypernode, it's a leaf <=> it has no incoming edges
    for ( size_t i = 0; i < N; ++i) {
        auto elist = H->incident(i);
        if ( elist.size() == 0 ) {
            leafHypernodes.push_back(i);
        }
    }

    typedef typename GT::vertex_descriptor NodeT;
    typedef unordered_set<NodeT> NodeSetT;

    // Is the node e contained in the set s?
    auto contains = [] ( const NodeSetT & s, NodeT e ) {
        return s.find(e) != s.end();
    };

    NodeSetT extantNodes;

    Google<int>::Set leafIds;
    leafIds.set_empty_key(-1);

    for ( auto l : T->getLeavesId() ) {
        leafIds.insert(l);
    }

    auto vp = boost::vertices(G);
    for ( auto it = vp.first; it != vp.second; ++it ) {
        auto v = *it;
        auto idx = G[v].idx;
        // found this node's id in the set of extant vertices
        if ( leafIds.find(idx) != leafIds.end() ) {
            extantNodes.insert(v);
        }
    }

    // Map from tree node ID to graph vertex ID
    unordered_map<int, NodeT> idToVertMap;
    for ( auto v = boost::vertices(G).first; v != boost::vertices(G).second; ++ v ) {
        idToVertMap[ G[*v].idx ] = *v;
    }
    //cerr << "ID TO VERTEX MAP SIZE = " << idToVertMap.size() << "\n";
    size_t nlost = 0;
    size_t nef = 0;
    double tweight = 0.0;
    // For every leaf hypernode
    for ( auto n : leafHypernodes ) {
        auto nd = H->vertex(n);

        auto nameU = T->getNodeName(nd.u());
        auto nameV = T->getNodeName(nd.v());

        // Check to see if u, v, or both have been lost
        auto endOfMap = idToVertMap.end();
        bool lostU = ( idToVertMap.find( nd.u() ) == endOfMap ) || ( ! contains(extantNodes, idToVertMap[nd.u()]) );
        bool lostV = ( idToVertMap.find( nd.v() ) == endOfMap ) || ( ! contains(extantNodes, idToVertMap[nd.v()]) );

        // The cost to / between lost nodes is always 0
        if ( lostU || lostV ) {
            auto lostCost = 0.5;
            vector<size_t> ev;
            Google<Derivation::flipT>::Set es;
            es.set_empty_key( make_tuple(-1, -1, ""));
            nlost += 1;
            slnDict[n] = { {0, Derivation(lostCost, n, ev, es)} };
        } else {
            // Otherwise, u and v both exist in the extant
            // network, so get the appropriate info
            auto u = idToVertMap[nd.u()];
            auto v = idToVertMap[nd.v()];
            auto f = nd.f();
            auto r = nd.r();

            if (u != v) {
                EdgeT fedge, redge;
                bool d_f, d_r;
                tie(fedge, d_f) = edge(u, v, G);
                double w_f = d_f ? G[fedge].weight : 0.0;
                tie(redge, d_r) = edge(v, u, G);
                double w_r = d_r ? G[redge].weight : 0.0;
                if ( undirected ) {
                    assert( w_f == w_r );
                }
                auto cost = (d_f == nd.f() and d_r == nd.r()) ? 1.0 : 0.0;
                Google<Derivation::flipT>::Set effectiveEdges;
                effectiveEdges.set_empty_key( make_tuple(-1, -1, ""));
                vector<size_t> ev;
                slnDict[n] = { {0, Derivation(cost, n, ev, effectiveEdges)} };

            } else {
                EdgeT e;
                bool hasSelfLoop;
                tie(e, hasSelfLoop) = edge(u, v, G);
                double w_l = hasSelfLoop ? G[e].weight : 0.0;

                auto cost = ( nd.f() == hasSelfLoop ) ? 1.0 : 0.0;
                Google<Derivation::flipT>::Set effectiveEdges;
                effectiveEdges.set_empty_key( make_tuple(-1, -1, ""));
                vector<size_t> ev;
                slnDict[n] = { {0, Derivation(cost, n, ev, effectiveEdges)} };

            } // ( u != v )
        } // ( lostU || lostV )
    } // loop over leaf hypernodes
}

template< typename GT >
void leafCostDict( unique_ptr<ForwardHypergraph> &H, TreePtrT &T, GT &G, bool directed, double cc, double dc, slnDictT &slnDict ) {
    /*
      Given the duplication tree T, the root vertex rv, the extant graph G and
      the constraints, fill in the values for the leaf nodes of the hypergraph
    */
    typedef typename boost::graph_traits< GT >::edge_descriptor EdgeT;
    boost::timer::auto_cpu_timer timer;
    auto undirected = !directed;
    auto isLost = [&]( int nid ) -> bool { return (T->getNodeName(nid)).find("LOST") != std::string::npos; };
    // Cost of going from the state encoded by a node to the state of a true graph edge

    auto none = make_tuple(false, false);
    auto fw = make_tuple(true, false);
    auto rev = make_tuple(false, true);
    auto both = make_tuple(true, true);

    auto costDict = getCostDict(cc, dc, directed);
    auto selfLoopCostDict = getSelfLoopCostDict(cc, dc, directed);
    auto costFunDict = getCostFunDict(cc, dc, directed);
    auto selfLoopCostFunDict = getSelfLoopCostFunDict(cc, dc, directed );

    auto N = H->order();
    auto M = H->size();

    //cerr << "total # of hypernodes = " << N << "\n";
    //cerr << "total # of hyperedges = " << M << "\n";

    // The list of all hypernodes with no descendants
    // We'll have either N^2 or (N^2)/2 leaf hypernodes (depending
    // on directedness)
    auto numExtantNodes = std::distance( vertices(G).first, vertices(G).second );
    auto numConn = ( numExtantNodes * numExtantNodes );
    if (undirected) {
        numConn /= 2;
    }
    vector<size_t> leafHypernodes;
    leafHypernodes.reserve(numConn);

    // For every hypernode, it's a leaf <=> it has no incoming edges
    for ( size_t i = 0; i < N; ++i) {
        auto elist = H->incident(i);
        if ( elist.size() == 0 ) {
            leafHypernodes.push_back(i);
        }
    }

    //cerr << "# of leaf hypernodes = " << leafHypernodes.size() << "\n";

    typedef typename GT::vertex_descriptor NodeT;
    typedef unordered_set<NodeT> NodeSetT;

    // Is the node e contained in the set s?
    auto contains = [] ( const NodeSetT & s, NodeT e ) {
        return s.find(e) != s.end();
    };

    NodeSetT extantNodes;

    Google<int>::Set leafIds;
    leafIds.set_empty_key(-1);

    for ( auto l : T->getLeavesId() ) {
        leafIds.insert(l);
    }

    auto vp = boost::vertices(G);
    for ( auto it = vp.first; it != vp.second; ++it ) {
        auto v = *it;
        auto idx = G[v].idx;
        // found this node's id in the set of extant vertices
        if ( leafIds.find(idx) != leafIds.end() ) {
            extantNodes.insert(v);
        }
    }

    // Map from tree node ID to graph vertex ID
    unordered_map<int, NodeT> idToVertMap;
    for ( auto v = boost::vertices(G).first; v != boost::vertices(G).second; ++ v ) {
        idToVertMap[ G[*v].idx ] = *v;
    }
    //cerr << "ID TO VERTEX MAP SIZE = " << idToVertMap.size() << "\n";
    size_t nlost = 0;
    size_t nef = 0;
    double tweight = 0.0;
    // For every leaf hypernode
    for ( auto n : leafHypernodes ) {
        auto nd = H->vertex(n);

        auto nameU = T->getNodeName(nd.u());
        auto nameV = T->getNodeName(nd.v());

        // Check to see if u, v, or both have been lost
        auto endOfMap = idToVertMap.end();
        //bool lostU = isLost(nd.u());//( idToVertMap.find( nd.u() ) == endOfMap ) || ( ! contains(extantNodes, idToVertMap[nd.u()]) );
        //bool lostV = isLost(nd.v());//( idToVertMap.find( nd.v() ) == endOfMap ) || ( ! contains(extantNodes, idToVertMap[nd.v()]) );
        bool lostU = ( idToVertMap.find( nd.u() ) == endOfMap ) || ( ! contains(extantNodes, idToVertMap[nd.u()]) );
        bool lostV = ( idToVertMap.find( nd.v() ) == endOfMap ) || ( ! contains(extantNodes, idToVertMap[nd.v()]) );

        // The cost to / between lost nodes is always 0
        if ( lostU || lostV ) {
            auto lostCost = 0.0;
            vector<size_t> ev;
            Google<Derivation::flipT>::Set es;
            es.set_empty_key( make_tuple(-1, -1, ""));
            nlost += 1;
            slnDict[n] = { {0, Derivation(lostCost, n, ev, es)} };
        } else {
            // Otherwise, u and v both exist in the extant
            // network, so get the appropriate info
            auto u = idToVertMap[nd.u()];
            auto v = idToVertMap[nd.v()];
            auto f = nd.f();
            auto r = nd.r();

            if (u != v) {
                EdgeT fedge, redge;
                bool d_f, d_r;
                tie(fedge, d_f) = edge(u, v, G);
                double w_f = d_f ? G[fedge].weight : 0.0;
                tie(redge, d_r) = edge(v, u, G);
                double w_r = d_r ? G[redge].weight : 0.0;
                if ( undirected ) {
                    assert( w_f == w_r );
                }

                tweight += (w_f + w_r) * (d_f + d_r) * (nd.f() + nd.r());

                //auto costFlipProb = costDict[ make_tuple(f,r) ][ make_tuple(d_f,d_r) ];
                auto costFlipProb = costFunDict[ make_tuple(f, r) ][ make_tuple(d_f, d_r) ](w_f, w_r);

                /*
                if (costFlip != costFlipProb) {
                    cerr << "whoops for transition (" << f << ", " << r << ") => (" << d_f << ", " << d_r << "), and (w_f, w_r) = (" << w_f << ", " << w_r << ")\n";
                    cerr << "costFlip = (" << get<0>(costFlip) << ", " << get<1>(costFlip) << "), but costFlipProb = (" << get<0>(costFlipProb) << ", " << get<1>(costFlipProb) << ")\n";
                    exit(1);
                }
                */

                auto cost = round3(get<0>(costFlipProb));
                auto flip = get<1>(costFlipProb);

                Google<Derivation::flipT>::Set effectiveEdges;
                effectiveEdges.set_empty_key( make_tuple(-1, -1, ""));

                if ( flip != "n" ) {
                    nef += 1;
                    effectiveEdges.insert( make_tuple(nd.u(), nd.v(), flip) );
                }
                vector<size_t> ev;
                slnDict[n] = { {0, Derivation(cost, n, ev, effectiveEdges)} };

            } else {
                EdgeT e;
                bool hasSelfLoop;
                tie(e, hasSelfLoop) = edge(u, v, G);
                double w_l = hasSelfLoop ? G[e].weight : 0.0;

                tweight += w_l * (hasSelfLoop) * (nd.f() + nd.r());
                //auto costFlipProb = selfLoopCostDict[ make_tuple(f,r) ][ hasSelfLoop ];
                auto costFlipProb = selfLoopCostFunDict[ make_tuple(f, r) ][ hasSelfLoop ]( w_l );
                /*
                if (costFlip != costFlipProb) {
                    cerr << "whoops for self loop transition (" << f << ", " << r << ") => (" << hasSelfLoop << "), and (w_l) = (" << w_l << ")\n";
                    cerr << "costFlip = (" << get<0>(costFlip) << ", " << get<1>(costFlip) << "), but costFlipProb = (" << get<0>(costFlipProb) << ", " << get<1>(costFlipProb) << ")\n";
                    exit(1);
                }
                */

                auto cost = round3(get<0>(costFlipProb));
                auto flip = get<1>(costFlipProb);

                Google<Derivation::flipT>::Set effectiveEdges;
                effectiveEdges.set_empty_key( make_tuple(-1, -1, ""));

                if ( flip != "n" ) {
                    nef += 1;
                    effectiveEdges.insert( make_tuple(nd.u(), nd.v(), flip) );
                }
                vector<size_t> ev;
                slnDict[n] = { {0, Derivation(cost, n, ev, effectiveEdges)} };

            } // ( u != v )
        } // ( lostU || lostV )
    } // loop over leaf hypernodes

    double sum = 0.0;
    //cerr << "slnDict size = " << slnDict.size() << "\n";
    for ( auto n : slnDict ) {
        sum += (n.second)[0].cost;
    }
    /*
    cerr << "SUM = " << sum << "\n";
    cerr << "NUM LOST = " << nlost << "\n";
    cerr << "NUM EFFECTIVE EDGES = " << nef << "\n";
    cerr << "TOTAL WEIGHT = " << tweight << "\n";
    */
    /*
    auto fname = "BASECASES.txt";
    std::fstream output( fname, std::fstream::out | std::fstream::trunc );
    for( auto kvit : slnDict ) {
        auto vert = H->vertex(kvit.first);
        if (! (isLost(vert.u()) || isLost(vert.v())) ) {
            output << std::fixed << std::setprecision(18) << T->getNodeName(vert.u()) << "\t" << T->getNodeName(vert.v()) << "\t" <<
                    (vert.f() ? "true" : "false") << "\t" << (vert.r() ? "true" : "false") << "\t" <<
                    (kvit.second)[0].cost << "\n";
        }
    }
    */
}

template <typename CostClassT>
tuple<double, cl_I> getCostCount( vector<vector<CostClassT>> &tkd,
                                  const vector<size_t> &bp,
                                  const size_t &eid,
                                  unique_ptr<ForwardHypergraph> &H ) {

    auto edge = H->edge(eid);
    double cost = edge.weight();
    cl_I count(1);
    size_t i = 0;
    for ( auto & tailNode : edge.tail() ) {
        // the cost class index
        size_t cci = bp[i];
        // We sum the costs and multiply the counts
        cost += tkd[tailNode][cci].cost();
        count *= tkd[tailNode][cci].total();
        i++;
    }
    return make_tuple(cost, count);
}

template <typename CostClassT>
double getCost( vector<vector<CostClassT> > &tkd,
                const vector<size_t> &bp,
                const size_t &eid,
                unique_ptr<ForwardHypergraph> &H ) {
    auto edge = H->edge(eid);
    double cost = edge.weight();
    size_t i = 0;
    for ( auto & tailNode : edge.tail() ) {
        // the cost class index
        size_t cci = bp[i];
        // We sum the costs and multiply the counts
        cost += tkd[tailNode][cci].cost();
        i++;
    }
    return round3(cost);
}

template <typename CostClassT>
cl_I getCount( vector< vector<CostClassT> > &tkd,
               const vector<size_t> &bp,
               const size_t &eid,
               unique_ptr<ForwardHypergraph> &H ) {
    auto edge = H->edge(eid);
    cl_I count(1);
    size_t i = 0;
    for ( auto & tailNode : edge.tail() ) {
        // the cost class index
        size_t cci = bp[i];
        // We sum the costs and multiply the counts
        count *= tkd[tailNode][cci].total();
        i++;
    }
    return count;
}

// 0 : 20 23 25
// 1 : 15 16 34
// 2 : 10 12 13
// rank, #
// (0, 3000)
//
vector< tuple<double, cl_I> > countEdgeSolutions(
    const double &ecost,
    const vector<size_t> &tailNodes,
    countDictT &countDict,
    const size_t &k,
    bool printMe,
    unique_ptr<ForwardHypergraph> &H,
    TreePtrT &t ) {

    auto vstr = [&]( const FlipKey & vert ) -> string {
        auto uname = t->getNodeName(vert.u());
        auto vname = t->getNodeName(vert.v());
        if (uname > vname) {
            auto tmp = vname;
            vname = uname;
            uname = tmp;
        }
        auto fstr = vert.f() ? "true" : "false";
        auto rstr = vert.r() ? "true" : "false";
        return uname + "\t" + vname + "\t" + fstr + "\t" + rstr;
    };

    // product pointers
    std::vector< size_t > elemSizes;
    elemSizes.reserve(tailNodes.size());
    double cost = ecost;
    for ( const auto & t : tailNodes ) {
        elemSizes.push_back( countDict[t].size() );
        cost += get<0>(countDict[t].front());
    }

    vector<dvsT> pq(1, make_tuple(cost, vector<size_t>(tailNodes.size(), 0)));
    QueueCmp<dvsT> ord;


    std::function< double( const vector<size_t>& ) > computeScore = [&] ( const vector<size_t> &inds ) -> double {
        size_t numNodes = tailNodes.size();
        double cost = ecost;
        for ( size_t i = 0; i < numNodes; ++i ) {
            cost += get<0>(countDict[ tailNodes[i] ][ inds[i] ]);
        }
        return cost;
    };

    typedef tuple<double, cl_I> ccT;
    vector< ccT > edgeSlns;
    double epsilon = 0.0;//5e-1;
    /*
    double minElement = round3(cost - 0.5*epsilon);
    vector< ccT > edgeSlns;
    for( size_t i = 1; i < k+2; ++ i) {
        edgeSlns.push_back( make_tuple(minElement+((i*epsilon)/2.0), cl_I(0)) );
    }

    auto computeScoreBin = [&]( double score ) -> size_t {
        return static_cast<size_t>(std::floor( std::abs(round3(score)-minElement) / epsilon ));
    };
    */

    size_t numSln = 0;
    bool finished = false;
    while ( !pq.empty() && numSln <= k ) {
        // Get the next best solution score from the top of the queue
        double cost;
        vector<size_t> inds;
        std::tie(cost, inds) = pq.front();
        std::pop_heap( pq.begin(), pq.end(), ord );
        pq.pop_back();

        // Compute the number of ways we can obtain this solution
        cl_I numSlns(1);
        for ( size_t i = 0; i < inds.size(); ++i ) {
            if (printMe) {
                auto vert = H->vertex(tailNodes[i]);
                cerr << vstr( vert ) << "\n";
                cerr << "i = " << i << ", score = " << get<0>(countDict[ tailNodes[i] ][ inds[i] ]) << ", count = " << get<1>(countDict[ tailNodes[i] ][ inds[i] ]) << "\n";
            }

            numSlns *= get<1>(countDict[ tailNodes[i] ][ inds[i] ]);
        }
        /*
        // put this solution into our # sln dict
        auto fp = std::find_if( edgeSlns.begin(), edgeSlns.end(), [=]( const ccT& cc ) ->bool { return (fabs(cost - get<0>(cc)) <= epsilon);  } );
        if (fp == edgeSlns.end()) { // we didn't find a solution  of this score
        if (printMe) {
        auto ess = (edgeSlns.size() == 0) ? "true" : "false";
        auto sd = ((edgeSlns.size() > 0) && (abs(cost - get<0>(edgeSlns.back())) >= epsilon )) ? "true" : "false";
        cerr << "edgeSlns.size() == 0 " <<  ess << "\n";
        if ( edgeSlns.size() > 0 ) {
        cerr << "cost = " << cost << ", cost(edgeSlns.back() = " << get<0>(edgeSlns.back()) << "\n";
        cerr << "cost - cost(edgeSlns.back()) = " << (cost - get<0>(edgeSlns.back())) << "\n";
        cerr << "abs(cost - cost(edgeSlns.back())) = " << std::abs(cost - get<0>(edgeSlns.back())) << "\n";
        cerr << "epsilon = " << epsilon << "\n";
        }
        cerr << "abs(cost - get<0>(edgeSlns.back())) > epsilon ) " <<  sd << "\n";
        }
        */
        if ( edgeSlns.size() == 0 || (fabs(cost - get<0>(edgeSlns.back())) >= epsilon) ) {
            if (printMe) {
                cerr << "IN IF: edgeSlns = ";
                for ( auto es : edgeSlns) {
                    cerr << get<1>(es) << ", ";
                }
                cerr << "\t new guy (" << cost << ", " << numSlns << ")\n";
            }
            edgeSlns.push_back( make_tuple( cost, numSlns ) );
        } else { // we found a solution of this score
            if (printMe) {
                cerr << "IN ELSE: edgeSlns = ";
                for ( auto es : edgeSlns) {
                    cerr << get<1>(es) << ", ";
                }
                cerr << "\t new guy (" << cost << ", " << numSlns << ")\n";
            }
            if ( cost < get<0>(edgeSlns.back()) ) {
                get<0>(edgeSlns.back()) = cost;
            }
            get<1>(edgeSlns.back()) += numSlns;
            /*if( cost < get<0>(*fp) ) {
                get<0>(*fp) = cost;
            }
            get<1>(*fp) = get<1>(*fp) + numSlns;
            */
        }
        /*
        size_t ind = computeScoreBin(cost);
        if (ind >= edgeSlns.size() ) {
            finished = true;
        } else {
            get<1>(edgeSlns[ind]) += numSlns;
            Utils::appendNext( cost, inds, elemSizes, pq, ord, computeScore );
        }
        */
        Utils::appendNext( cost, inds, elemSizes, pq, ord, computeScore );
        numSln = edgeSlns.size();
    }
    /*
    size_t i = 0;
    cl_I zero(0);
    while ( i <= k && get<1>(edgeSlns[i]) > zero ) {
        i += 1;
    }
    i = std::min(i,k);
    edgeSlns.resize(i);
    */
    if ( edgeSlns.size() == k + 1 ) {
        edgeSlns.resize(k);
    }
    return edgeSlns;
    /*
    std::sort(edgeSlns.begin(),edgeSlns.end(),
        [&]( const tuple<double, cl_I>& e1, const tuple<double, cl_I>& e2 ) -> bool {
            return get<0>(e1) < get<0>(e2);
        }
    );
    if(edgeSlns.size() == k+4) { for (size_t i = 0; i < 4; ++i) { edgeSlns.pop_back(); } }
    */
}

vector<double> computeAlphasDouble( const double &bscale, const vector<tuple<double, cl_I>> &slnVec, size_t k, const cl_I &total ) {
    if (slnVec.size() == 0) {
        return vector<double>();
    }
    vector<double> scores;
    scores.reserve(slnVec.size());
    double bestScore = get<0>(slnVec.front());
    double worstScore = get<0>(slnVec.back());
    if ( bestScore == worstScore && slnVec.size() > 2 ) {
        cerr << "bestScore (" << bestScore << ") == worstScore (" << worstScore << ")" << "\n";
        cerr << "=== slnVec ===\n";
        for (const auto & e : slnVec) {
            cerr << "score = " << get<0>(e) << ", count = " << get<1>(e) << "\n";
        }
        std::abort();
    }
    double diff = (worstScore == bestScore) ? 1.0 : worstScore - bestScore; // std::max(0.01, worstScore - bestScore);
    size_t N = slnVec.size();

    //double scale = (mpfr::log( J ) - mpfr::log( I )).toDouble() / diff;
    //double scale = 2.0 * estimateThermodynamicBeta( slnVec, bestScore ); // (6.9*k) / diff; // (1.25 * k) / diff;
    //double scale = 1.8;
    double beta = bscale * (k);
    double scale = bscale /  diff;//bscale / diff;
    //double scale = (0.5 * k) / diff;//(2.0 * N) / diff;//(1.5*k) / diff;
    //std::cerr << " **** beta = " << scale << " **** \n";

    //auto floatPrec = cln::float_format(40);
    //cl_F sum = cl_float(0.0,floatPrec);
    double sum(0.0);
    size_t i = 0;
    for (const auto & e : slnVec) {
        //scale = std::pow(bscale,i);
        double a = std::exp( -std::abs( (bestScore - get<0>(e)) * scale) );
        scores.push_back( a ); //cl_float(get<1>(e) * cl_float(a,floatPrec) ) );
        sum += scores.back();
        i += 1;
    }

    auto invSum = 1.0 / sum;
    vector<double> alphas;
    alphas.reserve(slnVec.size());
    for (const auto & s : scores) {
        alphas.push_back( s * invSum );
    }
    return alphas;
}

void viterbiCount( unique_ptr<ForwardHypergraph> &H, TreePtrT &t, TreeInfo &ti, double penalty, const vector<size_t> &order,
                   slnDictT &slnDict, countDictT &countDict, const size_t &k,
                   const string &outputName, const vector<FlipKey> &outputKeys, const double &beta ) {

    // Compute the *weighted* probability of each edge being in
    // the top k distinct scoring solutions

    cpplog::FileLogger log( "log.txt", true );
    auto vstr = [&]( const FlipKey & vert ) -> string {
        auto uname = t->getNodeName(vert.u());
        auto vname = t->getNodeName(vert.v());
        if (uname > vname) {
            auto tmp = vname;
            vname = uname;
            uname = tmp;
        }
        auto fstr = vert.f() ? "true" : "false";
        auto rstr = vert.r() ? "true" : "false";
        return uname + "\t" + vname + "\t" + fstr + "\t" + rstr;
    };

    double costsum = 0.0;
    // Each leaf has a single solution which is, by definition, of optimal cost
    for ( const auto & vit : order ) {
        if ( H->incident(vit).size() == 0 ) {
            countDict[vit].push_back( make_tuple(slnDict[vit][0].cost, cl_I(1)) );
            costsum += slnDict[vit][0].cost;
        }
    }


    auto cdname = "CDICT.txt";
    std::fstream cdout( cdname, std::fstream::out | std::fstream::trunc );
    cdout << countDict.size() << "\n";
    for ( auto ele : countDict ) {
        size_t vit = ele.first;
        auto vert = H->vertex(vit);
        cdout << vstr(vert) << "\t" << ele.second.size() << "\n";
        for ( auto tup : ele.second ) {
            cdout << std::fixed << std::setprecision(18) << get<0>(tup) << "\t" << get<1>(tup) << "\n";
        }
    }

    cerr << "COSTSUM = " << costsum << "\n";
    cerr << "ORDER SIZE = " << order.size() << "\n";
    cerr << "SIZE OF COUNT DICT = " << countDict.size() << "\n";
    typedef size_t edgeIdT;
    // For each edge, count the number of solutions having each score
    unordered_map< edgeIdT, unordered_map< double, cl_I > > edgeCountMap;
    unordered_map< edgeIdT, unordered_map< double, double > > edgeProbMap;

    // A map holding which edges are used to obtain *an* optimal
    // solution for each vertex
    unordered_map<size_t, unordered_map<double, unordered_set<size_t> >> usedEdges;
    //usedEdges.set_empty_key(std::numeric_limits<size_t>::max());

    auto N = order.size();//H->order();
    size_t ctr = 0;

    auto oname = "EDGES.txt";
    std::fstream out( oname, std::fstream::out | std::fstream::trunc );

    // For each vertex, in topological order (up)
    for ( auto vit = order.begin(); vit != order.end(); ++vit, ++ctr ) {
        if ( !(ctr % 1000) || ctr == N - 1 ) {
            cerr << "\r\rProcessed " << 100.0 * (static_cast<float>(ctr) / N) << "% of the vertices";
        }
        auto vert = H->vertex(*vit);

        // Put the results in an ordered map -- the sorted
        // property will be useful later for maintaining only the
        // k-best score classes
        map< double, vector<tuple<size_t, cl_I> > > edgeCostMap;

        //cerr << "SLN FOR NODE " << vstr(vert) << "\n";
        out << vstr(vert);// << H->incident(*vit).size() << "\n";
        std::vector<cl_I> nedgesln;
        // loop over all incoming edges and compute the # of
        // solutions over each edge as well as that solution's cost
        for ( const auto & e : H->incident(*vit) ) {

            auto edge = H->edge(e);
            auto hind = edge.head();
            auto w = edge.weight();
            //auto tvert = H->vertex(edge.tail()[0]);

            auto nameU = t->getNodeName(vert.u());
            auto nameV = t->getNodeName(vert.v());
            if ( nameU > nameV ) {
                std::swap(nameU, nameV);
            }
            bool printMe = false;//( nameU == "129064_Fr" );// && nameV == "n9758" );
            auto currentEdgeSlns = countEdgeSolutions( w, edge.tail(), countDict, k, printMe, H, t );
            if (printMe) {
                cerr << "EDGE: " << e << " counts = ";
                for ( auto el : currentEdgeSlns ) {
                    cerr << get<1>(el) << ", ";
                }
                cerr << "\n";
            }
            //if ( t->getNodeName(vert.u()) == "10866" )

            cl_I edgeSum(0);
            for ( const auto & ent : currentEdgeSlns ) {
                edgeSum += get<1>(ent);
            }
            nedgesln.push_back(edgeSum);

            out << w << "\t" << currentEdgeSlns.size();
            size_t j = 0;
            for ( const auto & ent : currentEdgeSlns ) {
                out << "\t" << get<0>(ent) << "\t" << get<1>(ent);
                j++;
            }
            out << "\n";

            /*
            if (ctr == N-1) {
                auto vert = H->vertex(*vit);
                auto un = t->getNodeName(vert.u());
                auto vn = t->getNodeName(vert.v());

                cout << "[" << un << ", " << vn << "] : (" << vert.f() << ", " << vert.r() << ")\n";
               for ( const auto& ent : currentEdgeSlns ) {
                 cout << "[score: " << get<0>(ent) << ", count: " << get<1>(ent) << "]\t";
               }
               cout << "\n";
            }
            */
            for ( const auto & ent : currentEdgeSlns ) {
                double score;
                cl_I count;
                tie(score, count) = ent;
                auto edgeContrib = make_tuple(e, count);
                edgeCostMap[score].push_back( edgeContrib );
                edgeCountMap[ e ][ score ] = count;
            }
        }

        std::sort(nedgesln.begin(), nedgesln.end(),
        []( const cl_I & x, const cl_I & y) {
            return x < y;
        }
                 );

        for ( auto e : nedgesln ) {
            out << "\t" << e;
        }
        out << "\n";

        // If we traversed any edges
        if ( edgeCostMap.size() > 0 ) {
            typedef tuple<double, cl_I, size_t> edgeSlnT;
            double minCost = std::numeric_limits<double>::max();
            size_t mk = std::min( edgeCostMap.size(), k );
            size_t ectr = 0;

            // for all incoming score classes
            for ( auto cmIt = edgeCostMap.begin(); cmIt != edgeCostMap.end() && ectr < mk; ++cmIt, ++ectr ) {
                // the score
                auto score = cmIt->first;
                // the set of edges providing this score
                const auto &providingEdges = cmIt->second;
                // the minimum cost incoming score
                minCost = std::min(minCost, score);
                // will count the number of solutions in this
                // score class
                cl_I numSln(0);

                // Update the information at the derived vertices
                for ( const auto & edgeCount : providingEdges ) {
                    size_t edgeInd;
                    cl_I count;
                    tie(edgeInd, count) = edgeCount;
                    usedEdges[*vit][score].insert( edgeInd );
                    // update the total number of solutions of the
                    // derived vertex
                    numSln += count;
                }
                // There are 'numSln' derivations yielding *vit at
                // a score of 'score'
                countDict[*vit].push_back( make_tuple(score, numSln) );

                // Now that we have a total count for the derived
                // vertex, compute each edge's contribution
                for ( const auto & edgeCount : providingEdges ) {
                    size_t edgeInd;
                    cl_I count;
                    tie(edgeInd, count) = edgeCount;
                    double edgeProb = double_approx( count / numSln );
                    // The probability of traversing this edge for
                    // derivations with the given score
                    edgeProbMap[edgeInd][score] = edgeProb;
                }

            }
            // Find the minimum cost edge
            Google<Derivation::flipT>::Set fs;
            fs.set_empty_key( Derivation::flipT(-1, -1, "") );
            slnDict[*vit][0] = Derivation( minCost, 0, vector<size_t>(), fs);

        }
    }
    // loop over verts
    cerr << "\n";

    auto epname = "EPROBS.txt";
    std::fstream epout( epname, std::fstream::out | std::fstream::trunc );
    for ( auto & em : edgeProbMap ) {
        //edgeProbMap.foreach{ case (eind, smap) =>
        auto eind = em.first;
        auto smap = em.second;
        epout << vstr( H->vertex(H->getHead(eind)) )  << "\t";
        epout << H->edge(eind).weight() << "\t" << smap.size()  << "\n";
        for ( auto sp : smap ) {
            epout << sp.first << "\t" << sp.second << "\n";
        }
    }

    typedef Google< size_t, double >::Map probMapT;

    auto getOrElse = [] ( probMapT & pm, const size_t & key, double alt ) {
        return (pm.find(key) == pm.end()) ? alt : pm[key];
    };

    cerr << "Backward step \n";

    probMapT probMap;
    probMapT outProbMap;
    probMap.set_empty_key( std::numeric_limits<size_t>::max() );
    outProbMap.set_empty_key( std::numeric_limits<size_t>::max() );

    FlipKey rootKey( t->getRootId(), t->getRootId(), false, false);
    auto rootInd = H->index(rootKey);
    // The root gets a probability of 1
    probMap[rootInd] = 1.0;
    //outProbMap[rootInd] = 1.0;

    ctr = 0;
    size_t tot = order.size();

    auto rootFlip = flipBoth(rootKey);
    H->addVertex( rootFlip );
    auto rootIdNoFlip = H->index(rootKey);
    auto rootIdFlip = H->index(rootFlip);

    // Compute the probabilities (outside)
    // Loop over all vertices in reverse topological order
    for ( auto vit = order.rbegin(); vit != order.rend(); ++vit, ++ctr ) {
        cerr << "\r\rprocessing node " << ctr << "/" << tot;
        // The current vertex and it's probability
        auto key = H->vertex(*vit);
        auto parentProb = 0.0;
        if ( probMap.find(*vit) != probMap.end() ) {
            parentProb = probMap[*vit];
        }

        // The total # of derivations of this vertex (over all
        // considered score classes)
        cl_I total(0);
        for ( size_t i = 0; i < countDict[*vit].size(); ++i ) {
            total += get<1>(countDict[*vit][i]);
        }

        // Compute the weights for each of the score classes
        // considered at this node
        auto alphas = computeAlphasDouble( beta, countDict[*vit], k, total );

        // for each top-k score class
        for ( size_t i = 0; i < countDict[*vit].size(); ++i ) {
            double pScore;
            cl_I pCount;
            // The score and it's count
            tie(pScore, pCount) = countDict[*vit][i];

            double tprob = 0.0;
            // for all incoming edges contributing to this score
            for ( const auto & e : usedEdges[*vit][pScore] ) {
                // The conditional probability of taking edge 'e'
                // given than we're deriving vertex *vit at score 'pScore'
                auto condProb = edgeProbMap[e][pScore];
                tprob += condProb;
                auto tail = H->getTail(e);

                auto ft = flipType( H->vertex(*vit), H->vertex(tail[0]) );
                auto outKey = keyForAction( H->vertex(*vit), ft );
                auto outInd = H->index(outKey);
                if ( outProbMap.find(outInd) == outProbMap.end() ) {
                    outProbMap[outInd] = 0.0;
                }
                outProbMap[outInd] += ( parentProb * (alphas[i] * condProb));

                // for all tail vertices of this edge
                for ( const auto & tind : tail ) {
                    // For each vertex in the tail of e, it gets
                    // probability mass for deriving *vit
                    // proportional to e's contribution
                    if ( probMap.find(tind) == probMap.end() ) {
                        probMap[tind] = 0.0;
                    }
                    probMap[tind] += (parentProb * ( alphas[i] * condProb ));
                    //outProbMap[tind] += (parentProb * ( alphas[i] * condProb ));
                }
            }

#ifdef DEBUG
            if ( usedEdges[*vit][pScore].size() > 0 && std::abs(tprob - 1.0) > 1e-4 ) {
                cerr << "ERROR : \n";
                cerr << "cond. probs from [" << t->getNodeName(key.u()) << ", " << t->getNodeName(key.v()) << "] (" << key.f() << ", " << key.r() << ")\n";
                cerr << "score = " << pScore << ", count = " << pCount << ", tprob = " << tprob << "\n";
                exit(1);
            }
#endif
        }
    }


    size_t i = 0;
    for ( const auto & sc : countDict[rootInd] ) {
        cout << "score class " << i << " has " << get<1>(sc) << " solutions of score  " << std::fixed << std::setprecision(16) << get<0>(sc) << "\n";
        ++i;
    }


    //bool restrictOutput = (outputKeys.size() != 0);
    bool restrictOutput = false;
    string fname = outputName;
    std::fstream output( fname, std::fstream::out | std::fstream::trunc );

    for ( size_t vid = 0; vid < H->order(); ++vid ) { // order.rbegin(); vit != order.rend(); ++vit ) {
        auto key = H->vertex(vid);
        if ( H->incident(vid).size() == 0 && vid != rootIdFlip && vid != rootIdNoFlip ) {
            if ( outProbMap[vid] != probMap[vid] && outProbMap[vid] > probMap[vid] ) {
                cout << "inProbMap has " << probMap[vid] << ", outProbMap has" << outProbMap[vid] << "\n";
            }
            outProbMap[vid] = probMap[vid];
        }

        auto approxInProb = probMap[vid];
        auto approxProb = outProbMap[vid];

        // ====================
        /*
        auto vert = H->vertex(*vit);
        std::string from, to;
        cout << "v = [" << t->getNodeName(vert.u()) << ", " << t->getNodeName(vert.v()) << "] (" << vert.f() << ", " << vert.r() << "), inProb = " << approxInProb << "\n";
        if ( !t->isLeaf(vert.u()) ) {
            auto rnode = vert.u();
            int LRN = t->getSonsId(rnode)[0]; int RRN = t->getSonsId(rnode)[1];
            cout << "children of " << t->getNodeName(vert.u()) << " are " << t->getNodeName(LRN) << " and " << t->getNodeName(RRN) << "\n";
        }

        if ( !t->isLeaf(vert.v()) ) {
            auto rnode = vert.v();
            int LRN = t->getSonsId(rnode)[0]; int RRN = t->getSonsId(rnode)[1];
            cout << "children of " << t->getNodeName(vert.v()) << " are " << t->getNodeName(LRN) << " and " << t->getNodeName(RRN) << "\n";

        }

        cout << "actions\n";
        auto parentProb = probMap[*vit];

        size_t tne = 0;
        // The total # of derivations of this vertex (over all
        // considered score classes)
        cl_I total(0);
        for ( size_t i = 0; i < countDict[*vit].size(); ++i ) { total += get<1>(countDict[*vit][i]); }

        // Compute the weights for each of the score classes
        // considered at this node
        auto alphas = computeAlphasDouble( countDict[*vit], k, total );


        // for each top-k score class
        for ( size_t i = 0; i < countDict[*vit].size(); ++i ) {
            double pScore; cl_I pCount;
            // The score and it's count
            tie(pScore, pCount) = countDict[*vit][i];

            cout << "score class " << i << "\n";

            double p = 0.0;
            // for all incoming edges contributing to this score
            for ( const auto& e : usedEdges[*vit][pScore] ) {
                // The conditional probability of taking edge 'e'
                // given than we're deriving vertex *vit at score 'pScore'
                auto condProb = edgeProbMap[e][pScore];
                auto prob = ( parentProb * (alphas[i] * condProb));
                auto tail = H->getTail(e);
                cout << "Action = " << flipType( key, H->vertex(tail[0])) << ", prob = " << prob << "\n";
                p += prob;
            }
            tne += usedEdges[*vit][pScore].size();
            cout << "score = " << pScore << ", numEdges =  " << usedEdges[*vit][pScore].size() << ", prob = " << p << "\n";
        }

        cout << "total number of incoming edges = " << H->incident(*vit).size() << ", total num used edges = " << tne << ", ";
        cout << "outProb = " << approxProb << "\n\n";
        */
        // ====================
        bool writeOut = true;
        if ( restrictOutput ) {
            writeOut = ( std::find( outputKeys.begin(), outputKeys.end(), key ) != outputKeys.end() );
        }

        if ( approxProb > 0.0 && writeOut ) {
            auto fs = flipStrMap[key.getDirTuple()];
            if ( fs != "n" ) {
                output << t->getNodeName(key.u()) << "\t" << t->getNodeName(key.v())
                       << "\t" << fs << "\t" << approxProb << "\n";
            }
        }
    }
    output.close();
    std::cerr << "\n";

}

double computeVertexProbability(const size_t &vid,
                                const TreePtrT& t,
                                const std::vector<double> &probs,
                                const std::vector<bool> &normed,
                                unique_ptr<ForwardHypergraph> &H,
                                Model &model) {

    auto printWithNames = [&] (size_t ind) -> void {
        auto k = H->vertex(ind);
        auto dirStr = (k.f() and k.r()) ? " <--> " : " X ";
        std::cerr << "[ " << t->getNodeName(k.u()) << ", " << t->getNodeName(k.v()) << 
            " : " << dirStr << "] ";
    };

    auto parentVertex = H->vertex(vid);
    //std::cerr << "\n\nP("; printWithNames(vid); std::cerr << ") = ";
    double prob = 0.0;
    size_t k = 0;
    // For each edge incident to this vertex
    for (auto & eid : H->incident(vid)) {
        //std::cerr << "[";
        // probability of the tail vertices
        auto edge = H->edge(eid);
        double tprob = 1.0;
        //std::cerr << "(";
        size_t i = 0;
        for (auto tid : edge.tail()) {
            auto childVertex = H->vertex(tid);
            //std::cerr << "P ("; printWithNames(tid); std::cerr << ") = ";
            auto cprob = probs[tid] * model.transitionProbability( childVertex, parentVertex );
            /*std::cerr << cprob;
            if ( i < edge.tail().size() - 1 ) {
                std::cerr << " * ";
                ++i;
            }
            assert( normed[tid] );*/
            tprob *= probs[tid] * model.transitionProbability( childVertex, parentVertex );
        }
        /*std::cerr << ")";
        std::cerr << "]"; if ( k != H->incident(vid).size() - 1 ) {
            std::cerr << " + ";
            ++k;
        }*/
        prob += tprob;
    }
    //std::cerr << " = " << prob << "\n\n";
    return prob;

}

template <typename CostClassT>
vector<CostClassT> computeKBest(const size_t &vid,
                                const size_t &k,
                                vector< vector<CostClassT> > &tkd,
                                unique_ptr<ForwardHypergraph> &H) {

    // Dictionary that holds, for each incoming edge, the
    // number of score classes for each tail node
    unordered_map<size_t, vector<size_t>> sizeDict;


    // Priority queue of derivations for the given vertex
    using boost::heap::fibonacci_heap;
    QueueCmp<edvsT> ord;
    fibonacci_heap<edvsT, boost::heap::compare<QueueCmp<edvsT>>> vpq(ord);

    // Will hold the top-k cost classes for this solution
    vector<CostClassT> cc;
    cc.reserve(k);

    for ( auto & eid : H->incident(vid) ) {
        // The edge, backpointer array and cost of the derivation
        auto edge = H->edge(eid);
        vector<size_t> bp(edge.tail().size(), 0);
        auto cost = round3(getCost(tkd, bp, eid, H));

        // Push this potential derivation on the queue
        vpq.push( make_tuple( cost, eid, bp ) );

        // Fill in the size array for this edge
        vector<size_t> sizes;
        sizes.reserve(edge.tail().size());
        for ( auto & tn : edge.tail() ) {
            sizes.push_back(tkd[tn].size());
        }
        sizeDict[eid] = sizes;
    }

    // Compute the cost of a derivation
    std::function< double( const size_t &, const vector<size_t>& ) > computeScore = [&] (const size_t & eid,
    const vector<size_t> &inds ) -> double {
        return round3( getCost(tkd, inds, eid, H) );
    };

    // Exact score classes
    double epsilon = 0.1;

    // While there are still derivations left in the heap,
    // and we don't yet have the required number of solutions
    while ( !vpq.empty() && cc.size() <= k ) {

        // Get the next best solution score from the top of the queue
        double cost;
        size_t eid;
        vector<size_t> inds;
        std::tie(cost, eid, inds) = vpq.top();
        vpq.pop();

        // Create the derivation
        auto count = getCount( tkd, inds, eid, H );
        CountedDerivation cderiv( cost, eid, inds, count );

        // If the list of cost classes is empty, or if this
        // derivation is belongs in a new cost class
        if ( cc.size() == 0 || (fabs(cost - cc.back().cost())) > epsilon ) {
            // Create the new cost class and append the
            // counted derivation
            CostClassT cclass(cost);
            cclass.appendToCount(cderiv);
            cc.push_back( cclass );
        } else { // we found a solution of this score
            // Append the counted derivation to the last cost class
            cc.back().appendToCount(cderiv);
        }

        // Append the successors of this derivation to the heap
        Utils::appendNextWithEdge( eid, inds, sizeDict[eid], vpq, computeScore );
    }

    // The cost class, the whole cost class and nothing
    // but the cost class
    if ( cc.size() == k + 1 ) {
        cc.resize(k);
    }
    return cc;
}


template <typename CostClassT>
void probabilistic( unique_ptr<ForwardHypergraph> &H,
                    Model &model,
                    TreePtrT &t,
                    const vector<size_t> &order,
                    slnDictT &slnDict,
                    const string &outputName,
                    const vector<FlipKey> &outputKeys ) {

    // Compute the *weighted* probability of each edge being in
    // the top k distinct scoring solutions
    cpplog::FileLogger log( "log.txt", true );

    size_t maxID = *(std::max_element(order.begin(), order.end()));

    // Dictionary that holds the top-k cost classes for each
    // vertex, as well as other relevant information
    vector< vector<CostClassT> > tkd(maxID + 1, vector<CostClassT>() ); //unordered_map< size_t, vector<CostClass> > tkd;

    /**
    * Dictionary that holds the probability of each state for each vertex
    */
    typedef double LogProbT;
    vector< LogProbT > probs( maxID + 1, 0.0 );
    vector< bool > normed( maxID + 1, false );

    auto vstr = [&]( const FlipKey & vert ) -> string {
        auto uname = t->getNodeName(vert.u());
        auto vname = t->getNodeName(vert.v());
        if (uname > vname) {
            auto tmp = vname;
            vname = uname;
            uname = tmp;
        }
        auto fstr = vert.f() ? "true" : "false";
        auto rstr = vert.r() ? "true" : "false";
        return uname + "\t" + vname + "\t" + fstr + "\t" + rstr;
    };


    auto printWithNames = [&] (size_t ind) -> void {
        auto k = H->vertex(ind);
        auto dirStr = (k.f() and k.r()) ? " <--> " : " X ";
        std::cerr << "[ " << t->getNodeName(k.u()) << ", " << t->getNodeName(k.v()) << 
            " : " << dirStr << "] ";
    };
    
    double costsum = 0.0;

    // Each leaf has a single solution which is, by definition, of optimal cost
    for ( const auto & vit : order ) {
        if ( H->isLeaf(vit) ) {
            probs[vit] = slnDict[vit][0].cost;
            normed[vit] = true;
        }
    }

    typedef size_t edgeIdT;
    auto N = order.size();
    size_t ctr = 0;
    ProgressDisplay showProgress(order.size());

    // For each vertex, in topological order (up)
    for ( auto vit = order.begin(); vit != order.end(); ++vit, ++ctr, ++showProgress ) {

        if ( H->isInternal(*vit) ) {

         for ( auto e : H->incident(*vit) ) {
            auto edge = H->edge(e);
            for ( auto v : edge.tail() ) {
               auto key = H->vertex(v);
               auto opKey = flipBoth(key);
               auto opInd = H->index(opKey);
               if ( !H->isLeaf(v) ) {
                 if ( normed[v] or normed[opInd] ) { 
                    if(normed[v] and !normed[opInd]) {
                        printWithNames(v);
                        std::cerr << "was normed but ";
                        printWithNames(opInd);
                        std::cerr << " was not\n";

                    } else if (!normed[v] and normed[opInd]) {
                        printWithNames(opInd);
                        std::cerr << "was normed but ";
                        printWithNames(v);
                        std::cerr << " was not\n";
                    } 
                    continue; 
                }

                 auto probSum = (probs[v] + probs[opInd]);
                 if ( ! (probSum > 0) ) {
                    std::cerr << "P("; printWithNames( v ); std::cerr << ") = " << probs[v] << "\n";
                    std::cerr << "P("; printWithNames( opInd); std::cerr << ") = " << probs[opInd] << "\n";
                    std::abort();
                 }
                 if ( probSum > 0 ) {
                 probs[v] /= probSum;
                 probs[opInd] /= probSum;
                 assert( probs[v] + probs[opInd] > 0.99 && probs[v] + probs[opInd] < 1.01 );
                 }
                 normed[v] = true; normed[opInd] = true;
               }
           }
        }

            auto vert = H->vertex(*vit);
            auto probVert = computeVertexProbability(*vit, t, probs, normed, H, model);
            probs[ *vit ] = probVert;
        }

    } // loop over verts

    typedef Google< size_t, double >::Map probMapT;
    size_t invalidIdx = std::numeric_limits<size_t>::max();
    probMapT probMap;
    probMapT outProbMap;
    probMap.set_empty_key( invalidIdx );
    outProbMap.set_empty_key( invalidIdx );

    // Map from a vertex to the maximum cost class that is
    // considered in deriving any solution that is actually used.
    FlipKey rootKey( t->getRootId(), t->getRootId(), false, false);
    auto rootInd = H->index(rootKey);

    ctr = 0;
    size_t tot = order.size();
    auto rootFlip = flipBoth(rootKey);
    H->addVertex( rootFlip );
    auto rootIdNoFlip = H->index(rootKey);
    auto rootIdFlip = H->index(rootFlip);

    cerr << "Down phase\n";
    showProgress.restart(order.size());

    std::vector<LogProbT> outProbs( probs.size(), 0.0 );

    // Compute the probabilities (down)
    // Loop over all vertices in *reverse* topological order
    for ( auto vit = order.rbegin(); vit != order.rend(); ++vit, ++ctr, ++showProgress ) {
        if( !normed[*vit] ) {
            auto k = H->vertex(*vit);
            std::cerr << "key ";
            printWithNames(*vit);
            std::cerr << " is not normed!!\n";
            //std::abort();
        }
        /*
        // for all incoming edges contributing to this cost class
        for ( const auto & e : H->incident(*vit) ) {
          // What is the action along edge 'e' (i.e. how is
          // the state of the interaction function different
          // between *vit and the tail nodes of e)
          auto tail = H->edge(e).tail();
          auto ft = flipType( H->vertex(*vit), H->vertex(tail[0]) );
          auto outKey = keyForAction( H->vertex(*vit), ft );
          auto outInd = H->index(outKey);
          outProbs[outInd] += probs[*vit];
        }
        */
    }

    //bool restrictOutput = (outputKeys.size() != 0);
    bool restrictOutput = false;
    string fname = outputName;
    std::fstream output( fname, std::fstream::out | std::fstream::trunc );

    for ( size_t vid = 0; vid < H->order(); ++vid ) { // order.rbegin(); vit != order.rend(); ++vit ) {
        auto key = H->vertex(vid);
        auto opKey = flipBoth(H->vertex(vid));
        auto opInd = H->index(opKey);
        auto approxProb = probs[vid];
        auto approxOtherProb = probs[ H->index(opKey) ];
        
        if ( (approxProb + approxOtherProb < 0.99 || approxProb + approxOtherProb > 1.01) and
             (t->isLeaf(key.u()) and t->isLeaf(key.v())) ) {
            std::cerr << "key = "; printWithNames(vid); std::cerr << ", P(key) = " << approxProb << "\n";
            std::cerr << "okey = "; printWithNames(opInd); std::cerr << ", P(okey) = " << approxOtherProb << "\n";
            //std::abort();
        }
        
        bool writeOut = true;
        if ( restrictOutput ) {
            writeOut = ( std::find( outputKeys.begin(), outputKeys.end(), key ) != outputKeys.end() );
        }

        if ( approxProb > 0.0 && writeOut ) {
            auto fs = flipStrMap[key.getDirTuple()];
            if ( fs != "n" ) {
                output << t->getNodeName(key.u()) << "\t" << t->getNodeName(key.v())
                       << "\t" << fs << "\t" << approxProb << "\n";
            }
        }
    }
    output.close();
    std::cerr << "\n";

}


template <typename CostClassT>
void viterbiCountNew( unique_ptr<ForwardHypergraph> &H, TreePtrT &t, TreeInfo &ti, double penalty, const vector<size_t> &order,
                      slnDictT &slnDict, countDictT &countDict, const size_t &k,
                      const string &outputName, const vector<FlipKey> &outputKeys, const double &beta ) {

    // Compute the *weighted* probability of each edge being in
    // the top k distinct scoring solutions
    cpplog::FileLogger log( "log.txt", true );

    size_t maxID = *(std::max_element(order.begin(), order.end()));

    // Dictionary that holds the top-k cost classes for each
    // vertex, as well as other relevant information
    vector< vector<CostClassT> > tkd(maxID + 1, vector<CostClassT>() ); //unordered_map< size_t, vector<CostClass> > tkd;

    auto vstr = [&]( const FlipKey & vert ) -> string {
        auto uname = t->getNodeName(vert.u());
        auto vname = t->getNodeName(vert.v());
        if (uname > vname) {
            auto tmp = vname;
            vname = uname;
            uname = tmp;
        }
        auto fstr = vert.f() ? "true" : "false";
        auto rstr = vert.r() ? "true" : "false";
        return uname + "\t" + vname + "\t" + fstr + "\t" + rstr;
    };

    double costsum = 0.0;

    // Each leaf has a single solution which is, by definition, of optimal cost
    for ( const auto & vit : order ) {
        if ( H->incident(vit).size() == 0 ) {
            tkd[vit] = { CostClassT(slnDict[vit][0].cost) };
            CountedDerivation cderiv( slnDict[vit][0].cost, std::numeric_limits<size_t>::max(), vector<size_t>(), cl_I(1) );
            tkd[vit].back().appendToCount( cderiv );
            countDict[vit].push_back( make_tuple(slnDict[vit][0].cost, cl_I(1)) );
            costsum += slnDict[vit][0].cost;
        }
    }

    /*
    cerr << "COSTSUM = " << costsum << "\n";
    cerr << "ORDER SIZE = " << order.size() << "\n";
    cerr << "SIZE OF COUNT DICT = " << countDict.size() << "\n";
    */

    typedef size_t edgeIdT;

    // For each edge, count the number of solutions having each score
    unordered_map< edgeIdT, unordered_map< double, cl_I > > edgeCountMap;
    unordered_map< edgeIdT, unordered_map< double, double > > edgeProbMap;

    auto N = order.size();
    size_t ctr = 0;
    ProgressDisplay showProgress(order.size());

    // For each vertex, in topological order (up)
    for ( auto vit = order.begin(); vit != order.end(); ++vit, ++ctr, ++showProgress ) {

        if ( H->incident(*vit).size() != 0 ) {
            auto vert = H->vertex(*vit);
            auto costClasses = computeKBest(*vit, k, tkd, H);
            tkd[*vit] = costClasses;
        }

        // Find the minimum cost edge
        Google<Derivation::flipT>::Set fs;
        fs.set_empty_key( Derivation::flipT(-1, -1, "") );
        slnDict[*vit][0] = Derivation( tkd[*vit].front().cost(), 0, vector<size_t>(), fs);

    } // loop over verts

    typedef Google< size_t, double >::Map probMapT;

    size_t invalidIdx = std::numeric_limits<size_t>::max();

    probMapT probMap;
    probMapT outProbMap;
    probMap.set_empty_key( invalidIdx );
    outProbMap.set_empty_key( invalidIdx );

    // Map from a vertex to the maximum cost class that is
    // considered in deriving any solution that is actually used.
    vector<size_t> maxCostClass(maxID + 1, 0);
    FlipKey rootKey( t->getRootId(), t->getRootId(), false, false);
    auto rootInd = H->index(rootKey);

    // We always consider all k classes for the root
    // There may be less than k classes if we have enumerated all
    // of them
    maxCostClass[rootInd] = std::min(k, tkd[rootInd].size());

    // The root gets a probability of 1
    probMap[rootInd] = 1.0;

    ctr = 0;
    size_t tot = order.size();
    auto rootFlip = flipBoth(rootKey);
    H->addVertex( rootFlip );
    auto rootIdNoFlip = H->index(rootKey);
    auto rootIdFlip = H->index(rootFlip);

    cerr << "Down phase\n";
    showProgress.restart(order.size());

    // Compute the probabilities (down)
    // Loop over all vertices in *reverse* topological order
    for ( auto vit = order.rbegin(); vit != order.rend(); ++vit, ++ctr, ++showProgress ) {

        // The current vertex and it's probability
        auto key = H->vertex(*vit);
        auto parentProb = (probMap.find(*vit) == probMap.end()) ? 0.0 : probMap[*vit];

        // The total # of derivations of this vertex (over all considered cost classes)
        cl_I total(0);
        vector< tuple<double, cl_I> > cd;

        auto maxDeriv = maxCostClass[*vit];
        //auto maxDeriv = tkd[*vit].size();

        for ( size_t i = 0; i < maxDeriv; ++i ) {
            total += tkd[*vit][i].total();
            for ( auto & e : tkd[*vit][i].usedEdges() ) {
                auto frontier = tkd[*vit][i].getEdgeFrontier(e);
                auto tail = H->edge(e).tail();

                for ( size_t j = 0; j < tail.size(); ++j) {
                    auto tn = tail[j];
                    // eager
                    maxCostClass[tn] = std::max( frontier[j] + 1, maxCostClass[tn] );

                    // lazy
                    // maxCostClass[tn] = std::max( frontier[j][0]+1, maxCostClass[tn] );
                }
            }
            cd.push_back( make_tuple(tkd[*vit][i].cost(), tkd[*vit][i].total()) );
        }

        // Compute the weights for each of the score classes
        // considered at this node
        auto alphas = computeAlphasDouble( beta, cd, maxCostClass[*vit], total );
        // for each top-k score class
        for ( size_t i = 0; i < maxDeriv; ++i ) {
            // The i-th cost class for this node
            auto cc = tkd[*vit][i];
            double tprob = 0.0;

            // for all incoming edges contributing to this cost class
            for ( const auto & e : cc.usedEdges() ) {

                // The conditional probability of taking edge 'e'
                // given than we're deriving vertex *vit at score
                // the given cost
                auto condProb = cc.edgeProb(e);
                tprob += condProb;
                auto tail = H->getTail(e);

                // What is the action along edge 'e' (i.e. how is
                // the state of the interaction function different
                // between *vit and the tail nodes of e)
                auto ft = flipType( H->vertex(*vit), H->vertex(tail[0]) );
                auto outKey = keyForAction( H->vertex(*vit), ft );
                auto outInd = H->index(outKey);

                // If we currently have not accumulated any
                // probability mass at outInd, then set it to 0
                if ( outProbMap.find(outInd) == outProbMap.end() ) {
                    outProbMap[outInd] = 0.0;
                }

                // Accumulate the probability due to the current
                // cost class at outInd
                outProbMap[outInd] += ( parentProb * (alphas[i] * condProb));

                // For each tail vertex of this hyperarc,
                // accumulate probability mass at this vertex
                for ( const auto & tind : tail ) {
                    // The tail vertex gets probability mass for
                    // deriving *vit proportional to e's contribution
                    if ( probMap.find(tind) == probMap.end() ) {
                        probMap[tind] = 0.0;
                    }

                    probMap[tind] += (parentProb * ( alphas[i] * condProb ));
                }
            }
        }
    }
    cerr << "\n";

    size_t i = 0;
    for ( const auto & cc : tkd[rootInd] ) {
        cout << "cost class " << i << " has " << cc.total() << " solutions of score  " << std::fixed << std::setprecision(16) << cc.cost() << "\n";
        ++i;
    }

    //bool restrictOutput = (outputKeys.size() != 0);
    bool restrictOutput = false;
    string fname = outputName;
    std::fstream output( fname, std::fstream::out | std::fstream::trunc );

    for ( size_t vid = 0; vid < H->order(); ++vid ) { // order.rbegin(); vit != order.rend(); ++vit ) {
        auto key = H->vertex(vid);
        if ( H->incident(vid).size() == 0 && vid != rootIdFlip && vid != rootIdNoFlip ) {
            if ( outProbMap[vid] != probMap[vid] && outProbMap[vid] > probMap[vid] ) {
                LOG_WARN(log) << "inProbMap has " << probMap[vid] << ", outProbMap has" << outProbMap[vid] << "\n";
            }
            outProbMap[vid] = probMap[vid];
        }

        auto approxInProb = probMap[vid];
        auto approxProb = outProbMap[vid];

        bool writeOut = true;
        if ( restrictOutput ) {
            writeOut = ( std::find( outputKeys.begin(), outputKeys.end(), key ) != outputKeys.end() );
        }

        if ( approxProb > 0.0 && writeOut ) {
            auto fs = flipStrMap[key.getDirTuple()];
            if ( fs != "n" ) {
                output << t->getNodeName(key.u()) << "\t" << t->getNodeName(key.v())
                       << "\t" << fs << "\t" << approxProb << "\n";
            }
        }
    }
    output.close();
    std::cerr << "\n";

}

vector<cl_RA> computeAlphas( const vector<tuple<double, cl_I>> &slnVec ) {
    vector<cl_RA> scores;
    cl_RA invSum(0);
    cl_RA one(1);
    for ( const auto & e : slnVec) {
        scores.push_back( one / static_cast<size_t>(get<0>(e) + 1.0) );
        invSum += scores.back();
    }
    //cl_RA invSum = 1 / sum;
    // cerr << "invSum = " << invSum << "\n";
    vector<cl_RA> alphas;
    alphas.reserve(slnVec.size());
    for ( const auto & s : scores) {
        alphas.push_back( s / invSum );
    }
    return alphas;
}


double estimateThermodynamicBeta( const vector<tuple<double, cl_I>> &slnVecIn, const double &emin ) {
    if (slnVecIn.size() == 1) {
        return 0.0;
    }
    vector< tuple<double, mpreal> > slnVec;
    slnVec.reserve(slnVecIn.size());
    std::stringstream ss;
    for ( const auto & e : slnVecIn ) {
        ss << get<1>(e);
        mpreal n = ss.str();
        ss.clear();
        slnVec.push_back( make_tuple( get<0>(e), n ) );
    }
    double ctr = 0.0;
    double beta = 0.0;
    size_t n = slnVec.size();
    double scale = 1.0 / ( n - 1 );
    for ( auto i = slnVec.begin(); (i + 1) != slnVec.end(); ++i ) {
        auto j = i + 1;
        auto lnI = mpfr::log( get<1>(*i) );
        auto eI = get<0>(*i) - emin;
        auto lnJ = mpfr::log( get<1>(*j) );
        auto eJ = get<0>(*j) - emin;
        double denom = (eJ - eI);

        if ( denom >= 1.0 ) {
            beta += scale * ( (lnJ - lnI).toDouble() ) / denom;
            ctr += 1.0;
        }
    }

    //beta /= ctr;
    //std::cerr << " ===== beta = " << beta << " ===== \n";
    return beta;
}

vector<double> getAlphas( const vector< tuple<double, cl_I> > &slnVec, const cl_I &numPaths ) {
    typedef tuple<double, cl_I> elemT;

    auto minmax = std::minmax_element( slnVec.begin(), slnVec.end(), [] (const elemT & e0, const elemT & e1) {
        return get<0>(e0) < get<0>(e1);
    } );
    double bestScore = get<0>(*(minmax.first));
    double worstScore = get<0>(*(minmax.second));
    double diff = std::max(1.0, worstScore - bestScore);
    double scale = 2.0 / diff;
    double totalWeight = 0.0;
    double alpha = 0.0;
    vector<double> invScores;
    invScores.reserve(slnVec.size());
    for ( const auto & e : slnVec) {
        invScores.push_back( (alpha * double_approx(get<1>(e) / numPaths)) + ((1.0 - alpha) * std::exp( bestScore - get<0>(e) * scale ))  );
        totalWeight += invScores.back();
    }

    vector<double> alphas;
    alphas.reserve(slnVec.size());
    for ( const auto & e : invScores ) {
        alphas.push_back( e / totalWeight );
    }
    return alphas;
}



void insideOutside( unique_ptr<ForwardHypergraph> &H, TreePtrT &t, TreeInfo &ti, double penalty, const vector<size_t> &order, slnDictT &slnDict, countDictT &countDict, const string &outputName ) {

    // We'll use vectors for these since the nodes and vectors have a well
    // defined ordering.
    vector< tuple<double, cl_I> > insideScores( H->order() );
    vector< tuple<double, cl_I> > edgeWeightMap( H->size() );

    // Each leaf has a single solution which is, by definition, of optimal cost
    for ( const auto & vit : order ) {
        if ( H->incident(vit).size() == 0 ) {
            insideScores[vit] = make_tuple(slnDict[vit][0].cost, cl_I(1));
        } else {
            insideScores[vit] = make_tuple(0.0, cl_I(0));
        }
    }

    for ( size_t i = 0; i < H->size(); ++i ) {
        edgeWeightMap[i] = make_tuple(0.0, cl_I(0));
    }

    size_t ctr = 0;
    size_t tot = H->order();
    // For each vertex, in topological order (inside)
    for ( auto vit = order.begin(); vit != order.end(); ++vit, ++ctr ) {
        cerr << "\r\rprocessing node " << ctr << "/" << tot;
        vector< tuple<double, cl_I> > scores;
        scores.reserve( H->incident(*vit).size() );
        if ( H->incident(*vit).size() > 0 ) { // For every non-leaf
            cl_I totalPaths(1);
            // loop over all incoming edges
            for ( const auto & e : H->incident(*vit) ) {

                auto tail = H->getTail(e);
                vector< tuple<double, cl_I> > tailScores;
                tailScores.reserve( tail.size() );

                cl_I numPaths(1);
                cl_I summedNumPaths(0);
                // Loop over all tail vertices
                for ( const auto & tind : tail ) {
                    tailScores.push_back( insideScores[tind] );
                    numPaths *= get<1>(insideScores[tind]);
                    summedNumPaths += get<1>(insideScores[tind]);
                    cerr << "tailNode # paths = " << get<1>(insideScores[tind]) << "\n";
                }
                cerr << "summedNumPaths = " << summedNumPaths << "\n";
                totalPaths += numPaths;
                double avgTailCost = 0.0;
                size_t tind = 0;
                vector<double> talphas(getAlphas(tailScores, summedNumPaths));
                for ( const auto & s : tailScores ) {
                    //cerr << "(score, frac) = " << get<0>(s) << ", (" << talphas[tind] << ")\n";
                    avgTailCost += get<0>(s);// * (1.0 / tailScores.size());//talphas[tind];
                    tind += 1;
                }
                //if (std::isnan(avgTailCost) || ctr > 200 ) { std::abort(); }
                scores.push_back( make_tuple( H->edge(e).weight() + avgTailCost , numPaths ) );
            }

            double avgHeadScore = 0.0;
            for ( const auto & s : scores ) {
                avgHeadScore += get<0>(s) * 1.0 / scores.size();//double_approx( get<1>(s) / totalPaths );
            }

            vector<double> contribWeights(getAlphas(scores, totalPaths));
            insideScores[*vit] = make_tuple(avgHeadScore, totalPaths);

            // once we know the score of the parent node, we can
            // compute the relative contribution from each incoming edge
            size_t eind = 0;
            for ( const auto & e : H->incident(*vit) ) {
                edgeWeightMap[ e ] = make_tuple(contribWeights[eind], get<1>(scores[eind]));
                ++eind;
            }
        }
    }

    typedef Google< size_t, double >::Map probMapT;

    auto getOrElse = [] ( probMapT & pm, const size_t & key, double alt ) {
        return (pm.find(key) == pm.end()) ? alt : pm[key];
    };

    probMapT probMap;
    probMap.set_empty_key( std::numeric_limits<size_t>::max() );

    FlipKey rootKey( t->getRootId(), t->getRootId(), false, false);
    auto rootInd = H->index(rootKey);
    // The root gets a probability of 1
    probMap[rootInd] = 1.0;

    ctr = 0;
    // Loop over all vertices in reverse topological order
    for ( auto vit = order.rbegin(); vit != order.rend(); ++vit, ++ctr ) {
        cerr << "\r\rprocessing node " << ctr << "/" << tot;

        auto key = H->vertex(*vit);
        auto parentProb = probMap[*vit];

        for ( const auto & e : H->incident(*vit) ) {
            auto condProb = get<0>(edgeWeightMap[e]);
            // for all tail vertices of this edge
            auto tail = H->getTail(e);
            for ( const auto & tind : tail ) {
                probMap[tind] += (parentProb * condProb);
            }
        }
    }

    string fname = outputName;
    std::fstream output( fname, std::fstream::out | std::fstream::trunc );

    for ( auto vit = order.rbegin(); vit != order.rend(); ++vit ) {
        auto key = H->vertex(*vit);
        auto approxProb = probMap[*vit];//double_approx(probMap[*vit]);
        auto fs = flipStrMap[key.getDirTuple()];
        if ( approxProb > 0.0 && fs != "n" ) {
            output << t->getNodeName(key.u()) << "\t" << t->getNodeName(key.v())
                   << "\t" << fs << "\t" << approxProb << "\n";

        }
    }
    output.close();
    std::cerr << "\n";



}

FlipKey keyForAction( const FlipKey &fk , const string &ft ) {
    if ( ft == "n" ) {
        return fk;
    }
    if ( ft == "b+" || ft == "b-" ) {
        return flipBoth(fk);
    }
    if ( ft == "f+" || ft == "f-" ) {
        return flipForward(fk);
    }
    if ( ft == "r+" || ft == "r-" ) {
        return flipReverse(fk);
    }
}





void viterbi( unique_ptr<ForwardHypergraph> &H, TreePtrT &t, TreeInfo &ti, double penalty, const vector<size_t> &order, slnDictT &slnDict ) {
    auto N = H->order();
    size_t ctr = 0;
    for ( auto vit = order.begin(); vit != order.end(); ++vit, ++ctr ) {
        cerr << "\r\rProcessed " << 100.0 * (static_cast<float>(ctr) / N) << "% of the vertices";
        auto vert = H->vertex(*vit);
        //cout << "Processing : [" << t->getNodeName(vert.u()) << ", " << t->getNodeName(vert.v()) << " : (" << vert.f() << vert.r() << ")]\t";
        double cval = std::numeric_limits<double>::max();
        if ( slnDict.find(*vit) != slnDict.end() ) {
            cval = slnDict[*vit][0].cost;
        }

        //    cout << "COST = " << cval << "\n";

        for ( auto & e : H->incident(*vit) ) {
            auto edge = H->edge(e);
            auto hind = edge.head();
            auto w = edge.weight();
            auto tvert = H->vertex(edge.tail()[0]);

            Google<Derivation::flipT>::Set flips;
            flips.set_empty_key( make_tuple(-1, -1, ""));
            //set<Derivation::flipT> flips;
            string children;
            for ( auto tn : edge.tail() ) {
                auto cvert = H->vertex(tn);
                w += slnDict[tn][0].cost;
                /*std::set_union( slnDict[tn][0].flips.begin(), slnDict[tn][0].flips.end(),
                                slnDict[tn][0].flips.begin(), slnDict[tn][0].flips.end(),
                                std::inserter( flips, flips.begin() ) );
                */
                for ( auto & f : slnDict[tn][0].flips ) {
                    flips.insert(f);
                }

                //children += "[ " + t->getNodeName(cvert.u()) + ", " + t->getNodeName(cvert.v()) + " : (" + lexical_cast<string>(cvert.f()) + lexical_cast<string>(cvert.r()) +")] ";
            }

            auto ft = flipType(vert, tvert);
            if ( ft != "n" ) {
                flips.insert( make_tuple(vert.u(), vert.v(), ft ) );
            }

            if ( w <= cval) {
                vector<size_t> bp( edge.tail().size(), 0 );
                slnDict[hind][0] = Derivation(w, e, bp, flips);
                //cout << " updated cost from " << cval << " -> " << tval << " " << ((w > 0) ? "by flipping and" : "") <<  " using " << children << "\n";
                cval = w;
            }
        }
    }
    cerr << "\n";
}




/** Alg 3 from the paper */
//typedef std::pair<size_t, Derivation> TaggedDerivT;
//Google< tuple<size_t,size_t> >::Set recursedStore;
//Google<size_t, vector<Derivation>*>::Map cand;

CandStoreT cand;
DerivStoreT derivs;

/**
* Perform the initial pass over the hypergraph to compute the optimal score class
* for each of the vertices.
*/
std::vector<size_t> viterbiPass(
    unique_ptr<ForwardHypergraph> &H,
    DerivStoreT &derivs,
    vector<size_t> &order ) {

    std::vector<size_t> slnEdges( *std::max_element(order.begin(), order.end()), std::numeric_limits<size_t>::max() );
    for ( auto vit : order ) {
        // Visit the vertices in topological order

        // For each incoming edge
        for ( auto eit : H->incident(vit) ) {

            // Edge info
            auto edge = H->edge(eit);
            auto tailNodes = edge.tail(); // the tail nodes
            vector<size_t> bp(tailNodes.size(), 0); // backpointer
            // array (all 0's)
            double cost; cl_I count;
            std::tie (cost, count) = getCostCount(derivs, bp, eit, H);
            auto deriv = CountedDerivation(cost, eit, bp, count);

            /*if ( derivs.find(vit) == derivs.end() ) {
                derivs[vit] = vector<CostClass>();
                }*/

            double costClassWidth = 0.25;

            bool found = false;
            for ( auto & d : derivs[vit] ) {
                if ( fabs(d.cost() - cost) <= costClassWidth ) {
                    found = true;
                    d.appendToCount(deriv);
                }
            }
            if (!found) {
                LazyCostClass cc(cost); cc.appendToCount(deriv);
                derivs[vit].push_back(cc);
            }

        }
        // Throw away all but the top cost class for each vertex
        std::sort( derivs[vit].begin(), derivs[vit].end(),
        [](const LazyCostClass & a, const LazyCostClass & b) -> bool {
            return a.cost() < b.cost();
        });

        derivs[vit].resize(1);
        // Pick a random optimal edge
        if ( derivs[vit].back().usedEdges().size() > 0 ) {
            auto eid = derivs[vit].back().usedEdges().back();
            auto e = H->edge(eid);
            slnEdges[vit] = eid;
        }
    }
    return slnEdges;
}



//unordered_map< size_t, vector<Derivation>* > cand;

DerivStoreT &initKBest( unique_ptr<ForwardHypergraph> &H, std::vector<size_t> &order, slnDictT &slnDict ) {
    size_t maxID = *(std::max_element(order.begin(), order.end()));
    // Dictionary that holds the top-k cost classes for each
    // vertex, as well as other relevant information
    derivs = vector< vector<LazyCostClass> >(maxID + 1, vector<LazyCostClass>() );

    double costsum = 0.0;

    // Each leaf has a single solution which is, by definition, of optimal cost
    for ( const auto & vit : order ) {
        if ( H->incident(vit).size() == 0 ) {
            derivs[vit] = { LazyCostClass(slnDict[vit][0].cost) };
            CountedDerivation cderiv( slnDict[vit][0].cost, std::numeric_limits<size_t>::max(), vector<size_t>(), cl_I(1) );
            derivs[vit].back().appendToCount( cderiv );
            costsum += slnDict[vit][0].cost;
        }
    }

    return derivs;
}

template<typename T>
void printVector( const vector<T> &v ) {
    cerr << "[";
    for ( auto & e : v ) {
        cerr << e << " ";
    }
    cerr << "]";
}

void getCandidates(  unique_ptr<ForwardHypergraph> &H, size_t vid, size_t k, DerivStoreT &derivs ) {

    using boost::heap::fibonacci_heap;
    CountedDerivCmp<CountedDerivation> ord;
    cand[vid] = fibonacci_heap<CountedDerivation, boost::heap::compare<CountedDerivCmp<CountedDerivation>>>(ord);

    for ( auto e : H->incident(vid) ) {
        if (! derivs[vid].back().hasEdge(e) ) {
            auto edge = H->edge(e);
            auto tailNodes = edge.tail();
            auto tailSize = tailNodes.size();
            std::vector<size_t> bp(tailSize, 0);
            double cost;
            cl_I count;
            std::tie( cost, count ) = getCostCount(derivs, bp, e, H);
            auto deriv = CountedDerivation(cost, e, bp, count);
            cand[vid].push(deriv);
        }
    }
}

bool lazyNext(
    unique_ptr<ForwardHypergraph> &H,
    boost::heap::fibonacci_heap<CountedDerivation, boost::heap::compare<CountedDerivCmp<CountedDerivation>>> &localCandidates,
    size_t eind,
    const vector<size_t> &j,
    size_t kp,
    DerivStoreT &derivs) {

    auto edge = H->edge(eind);
    auto tailNodes = edge.tail();
    auto headNode = H->vertex(edge.head());
    auto edgeCost = edge.weight();

    size_t i = 0;
    while ( i < tailNodes.size() ) {
        vector<size_t> jp(j);
        jp[i] += 1;
        lazyKthBest(H, tailNodes[i], jp[i] + 1, kp, derivs);

        if ( jp[i] < derivs[tailNodes[i]].size() ) {
            double cost; cl_I count;
            std::tie(cost, count) = getCostCount(derivs, jp, eind, H);
            auto deriv = CountedDerivation(cost, eind, jp, count);
            localCandidates.push(deriv);
        }
        if (j[i] != 0) {
            return true;
        }
        i += 1;
    }
    return true;
}

bool lazyKthBest(
    unique_ptr<ForwardHypergraph> &H,
    size_t vid,
    size_t k,
    size_t kp,
    DerivStoreT &derivs) {


    // If this is our first vist to vertex v, then
    // populate its candidate list

    //cerr << "Searching for candidates for " << v << "\n";
    if ( cand.find(vid) == cand.end() ) {
        getCandidates(H, vid, kp, derivs);
    }

    // Until we have the required number of derivations
    while ( derivs[vid].size() <= k ) {
        // If we've already computed the last derivation

        if ( derivs[vid].size() > 0 ) {
            auto usedEdges =  derivs[vid].back().usedEdges();
            for ( auto e : usedEdges ) {
                auto &bps = derivs[vid].back().getEdgeFrontier(e);
                for ( auto & bp : bps ) {
                    lazyNext(H, cand[vid], e, bp, kp, derivs);
                }
            }
            // At this point, we no longer need the actual
            // derivations for this cost class, so free them
            derivs[vid].back().freeDerivations();
        }

        double costClassWidth = 1.0;//0.25;
        // Get the next best candidate from the heap
        if ( !cand[vid].empty() ) {
            auto topCost = cand[vid].top().cost;

            derivs[vid].push_back( LazyCostClass(topCost) );
            while ( !cand[vid].empty() && fabs(cand[vid].top().cost - topCost) <= costClassWidth ) {
                auto d = cand[vid].top();
                cand[vid].pop();
                derivs[vid].back().appendToCount(d);
            }

        } else {
            return false;
        }
    }

    return true;
}



void computePosteriors(
    unique_ptr<ForwardHypergraph> &H,
    const TreePtrT &t,
    vector<size_t> &order,
    DerivStoreT &derivs,
    const string &outputName,
    const vector<FlipKey> &outputKeys,
    const double &beta
) {


    typedef vector<InOutProb> ProbMapT;
    cerr << "Down step";
    ProbMapT probMap(H->order() + 2, InOutProb());


    cerr << "Size of probMap = " << probMap.size() << "\n";

    auto rootKey = FlipKey( t->getRootId(), t->getRootId(), false, false );
    auto rootInd = H->index(rootKey);

    // The root gets a probability of 1
    probMap[rootInd].inProb = 1.0;

    size_t ctr = 0;
    size_t tot = order.size();

    auto rootFlip = flipBoth(rootKey);
    H->addVertex(rootFlip);
    auto rootIdNoFlip = H->index(rootKey);
    auto rootIdFlip = H->index(rootFlip);

    for ( auto vitp = order.rbegin(); vitp != order.rend(); ++vitp ) {
        auto vit = *vitp;
        cerr << "\r\rprocessing node " << ctr << "/" << tot;
        ctr += 1;
        // The current vertex and it's probability
        auto key = H->vertex(vit);
        auto parentProb = probMap[vit].inProb;

        vector< tuple<double, cl_I> > cd;
        // The total # of derivations of this vertex (over all
        // considered score classes)
        cl_I total(0);
        for ( auto & cc : derivs[vit] ) {
            cd.push_back( make_tuple(cc.cost(), cc.total()) );
            total += cc.total();
        }

        // Compute the weights for each of the score classes
        // considered at this node
        auto alphas = computeAlphasDouble( beta, cd, derivs[vit].size(), total );

        // for each top-k score class
        for ( size_t i = 0; i < derivs[vit].size(); ++i ) {

            // The score and it's count
            double pScore(derivs[vit][i].cost());
            cl_I pCount(derivs[vit][i].total());

            double tprob = 0.0;
            // for all incoming edges contributing to this score
            auto incoming = derivs[vit][i].usedEdges();

            for ( auto e : incoming ) {


                // The conditional probability of taking edge 'e'
                // given than we're deriving vertex *vit at score 'pScore'
                double condProb = derivs[vit][i].edgeProb(e);
                tprob += condProb;
                auto tail = H->getTail(e);
                auto ft = flipType(H->vertex(vit), H->vertex(tail[0]));
                auto outKey = keyForAction(H->vertex(vit), ft);
                auto outInd = H->index(outKey);
                probMap[outInd].outProb += (parentProb * (alphas[i] * condProb));

                // for all tail vertices of this edge
                for ( auto tind : tail ) {

                    // For each vertex in the tail of e, it gets
                    // probability mass for deriving *vit
                    // proportional to e's contribution
                    probMap[tind].inProb += (parentProb * (alphas[i] * condProb));
                }
            }
        }
    }


    size_t i = 0;
    for ( const auto & cc : derivs[rootInd] ) {
        cout << "cost class " << i << " has " << cc.total() << " solutions of score  " << std::fixed << std::setprecision(16) << cc.cost() << "\n";
        ++i;
    }


    //bool restrictOutput = (outputKeys.size() != 0);

    bool restrictOutput = false;
    string fname = outputName;
    std::fstream output( fname, std::fstream::out | std::fstream::trunc );

    for ( size_t vid = 0; vid < H->order(); ++vid ) { // order.rbegin(); vit != order.rend(); ++vit ) {
        auto key = H->vertex(vid);
        if ( H->incident(vid).size() == 0 && vid != rootIdFlip && vid != rootIdNoFlip ) {
            probMap[vid].outProb = probMap[vid].inProb;
        }

        auto approxInProb = probMap[vid].inProb;
        auto approxProb = probMap[vid].outProb;

        bool writeOut = true;
        if ( restrictOutput ) {
            writeOut = ( std::find( outputKeys.begin(), outputKeys.end(), key ) != outputKeys.end() );
        }

        if ( approxProb > 0.0 && writeOut ) {
            auto fs = flipStrMap[key.getDirTuple()];
            if ( fs != "n" ) {
                output << t->getNodeName(key.u()) << "\t" << t->getNodeName(key.v())
                       << "\t" << fs << "\t" << approxProb << "\n";
            }
        }
    }
    output.close();
    std::cerr << "\n";
}

bool sameEffect( const Derivation &d0, const Derivation &d1 ) {
    if ( d0.cost != d1.cost ) {
        return false;
    }
    // otherwise, same scores, check the flips
    // vector<Derivation::flipT> res;
    // std::set_symmetric_difference( d0.flips.begin(), d0.flips.end(), d1.flips.begin(), d1.flips.end(), std::back_inserter(res) );
    // return (res.size() == 0) ? true : false;

    // for ( auto f0 : d0.flips ) {
    //   if ( d1.flips.find(f0) == d1.flips.end() ) { return false; }
    // }

    bool hasDiffElement = any_of( d0.flips.begin(),
                                  d0.flips.end(),
    [&] ( const tuple<int, int, string> &d ) {
        return d1.flips.find(d) == d1.flips.end();
    } );

    // for ( auto f1 : d1.flips ) {
    //     if ( d0.flips.find(f1) == d0.flips.end() ) { return false; }
    // }
    //cerr << "\n\nSame Effect:\n" << d0 << "\n" << d1 << "\n\n";
    return !hasDiffElement;
}






}

using GraphUtils::undirectedGraphT;
using GraphUtils::directedGraphT;
using std::unique_ptr;

template void MultiOpt::MLLeafCostDict< undirectedGraphT  >( unique_ptr<ForwardHypergraph> & , Utils::Trees::TreePtrT & , undirectedGraphT & , bool , double , double , MultiOpt::slnDictT & );

template void MultiOpt::MLLeafCostDict< directedGraphT >( unique_ptr<ForwardHypergraph> & , Utils::Trees::TreePtrT & , directedGraphT & , bool , double , double , MultiOpt::slnDictT & );

template void MultiOpt::leafCostDict< undirectedGraphT  >( unique_ptr<ForwardHypergraph> & , Utils::Trees::TreePtrT & , undirectedGraphT & , bool , double , double , MultiOpt::slnDictT & );

template void MultiOpt::leafCostDict< directedGraphT >( unique_ptr<ForwardHypergraph> & , Utils::Trees::TreePtrT & , directedGraphT & , bool , double , double , MultiOpt::slnDictT & );

using MultiOpt::LazyCostClass;
using MultiOpt::EagerCostClass;

template tuple<double, cl_I> MultiOpt::getCostCount( vector<vector<LazyCostClass>> &, const vector<size_t> &, const size_t &, unique_ptr<ForwardHypergraph> &);
template tuple<double, cl_I> MultiOpt::getCostCount( vector<vector<EagerCostClass>> &, const vector<size_t> &, const size_t &, unique_ptr<ForwardHypergraph> &);

template double MultiOpt::getCost( vector<vector<LazyCostClass>> &, const vector<size_t> &, const size_t &, unique_ptr<ForwardHypergraph> &);
template double MultiOpt::getCost( vector<vector<EagerCostClass>> &, const vector<size_t> &, const size_t &, unique_ptr<ForwardHypergraph> &);

template cl_I MultiOpt::getCount( vector<vector<LazyCostClass>> &, const vector<size_t> &, const size_t &, unique_ptr<ForwardHypergraph> &);
template cl_I MultiOpt::getCount( vector<vector<EagerCostClass>> &, const vector<size_t> &, const size_t &, unique_ptr<ForwardHypergraph> &);

template void MultiOpt::viterbiCountNew<CostClass<EdgeDerivInfoEager>>( unique_ptr<ForwardHypergraph> &H,
        Utils::Trees::TreePtrT &t,
        Utils::TreeInfo &ti,
        double penalty,
        const vector<size_t> &order,
        MultiOpt::slnDictT &slnDict,
        MultiOpt::countDictT &countDict,
        const size_t &k,
        const string &outputName,
        const vector<FlipKey> &outputKeys,
        const double &beta );

template void MultiOpt::probabilistic<CostClass<EdgeDerivInfoEager>>(
    unique_ptr<ForwardHypergraph> &H,
    Model &model,
    TreePtrT &t,
    const vector<size_t> &order,
    slnDictT &slnDict,
    const string &outputName,
    const vector<FlipKey> &outputKeys );
