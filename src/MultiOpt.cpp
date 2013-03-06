#include "MultiOpt.hpp"
#include "TreeUtils.hpp"
#include "ProgressDisplay.hpp"
#include "ezETAProgressBar.hpp"
#include "model.hpp"
//#include "mpl_cartprod.hpp"
#include "combination.hpp"

#include <stack>
#include <cmath>
#include <utility>
#include <cln/float.h>
#include <boost/timer/timer.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <boost/heap/pairing_heap.hpp>
#include <boost/heap/binomial_heap.hpp>
#include <boost/heap/skew_heap.hpp>
#include <boost/range/irange.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/copy.hpp>
#include <boost/range/counting_range.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <boost/pool/pool_alloc.hpp>
 
#include <Bpp/Phyl/TreeTools.h>

/** Google's dense hash set and hash map **/
#include <google/dense_hash_set>
#include <google/dense_hash_map>
#include <google/sparse_hash_map>

// For logging
#include "cpplog.hpp"

#ifdef DNDEBUG
#undef DNDEBUG
#include <cassert>
#define DNDEBUG
#endif

namespace MultiOpt {
using Utils::round3;

const auto none = make_tuple(false, false); //!< No edge exists between vertices u & v.
const auto fw = make_tuple(true, false);    //!< An edge (u,v) exists.
const auto rev = make_tuple(false, true);   //!< An edge (v,u) exists.
const auto both = make_tuple(true, true);   //!< If the graph is directed, then edges (u,v) & (v,u) exist;
                                            //!< otherwise, the edge {u,v} exists.

static const string flipStrMap[] = {"n", "f", "r", "b"};/*<! Map from a type of edge state to its textual representation */

/** \brief flipDict[A][B] gives the textual representation of action A -> B.
* Given edge states A & B, returns a textual representation of the action required to move from A to B (i.e. the action A->B). 
*/
static const string flipDict[][4] = 
{ 
  { "n", "f+", "r+", "b+" }, /*none transitions*/
  { "f-", "n", "f-r+", "r+" }, /*forward transitions*/
  { "r-", "f+r-", "n", "f+" }, /*reverse transitions*/
  { "b-", "r-", "f-", "n" }, /*both transitions*/
};

FlipState directionsToFlipState( bool f, bool r ) {
    if ( f and r ) { return FlipState::both; }
    if ( f and !r ) { return FlipState::forward; }
    if ( !f and r ) { return FlipState::reverse; }
    if ( !f and !r ) { return FlipState::none; }
}

costMapT getCostDict ( double cc, double dc, bool directed ) {
    auto undirected = ! directed;
    auto none = FlipState::none;
    auto fw = FlipState::forward;
    auto rev = FlipState::reverse;
    auto both = FlipState::both;

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
    auto none = FlipState::none;
    auto fw = FlipState::forward;
    auto rev = FlipState::reverse;
    auto both = FlipState::both;

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
    auto none = FlipState::none;
    auto fw = FlipState::forward;
    auto rev = FlipState::reverse;
    auto both = FlipState::both;

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
    auto none = FlipState::none;
    auto fw = FlipState::forward;
    auto rev = FlipState::reverse;
    auto both = FlipState::both;

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

/**
* Returns the textual representation of the action hvert -> tvert.
*/
inline string flipType( const FlipKey &hvert, const FlipKey &tvert ) {
    return flipDict[hvert.state()][tvert.state()];
}

template<typename GT>
std::unordered_set<int> projectToReversedGraph( unique_ptr<ForwardHypergraph> &H, GT &G ) {
    auto M = H->size();
    std::unordered_set<int> vset;

    size_t start = 0;
    for ( auto eid : boost::irange(start, M) ) { 
        auto e = H->edge(eid);
        auto h = e.head();
        vset.insert(h);
        for ( auto t : e.tail() ) {
            vset.insert(t);
            add_edge(h, t, G);
        }
    }

    return vset;
}

void topologicalOrder( unique_ptr<ForwardHypergraph> &H, TreePtrT& tree, const TreeInfo& ti, vector<size_t> &order ) {
    using boost::adjacency_list;
    using boost::vecS;
    using boost::directedS;
    using boost::undirectedS;

    cpplog::FileLogger log( "log.txt", true );
    auto isLost = [&]( int nid ) -> bool { return (tree->getNodeName(nid)).find("LOST") != std::string::npos; };
    typedef adjacency_list<vecS, vecS, directedS> graphT;
    typedef adjacency_list<vecS, vecS, undirectedS> undirGraphT;
    graphT G( H->order() );

    vector<size_t> torder;
    std::unordered_set<int> vset = projectToReversedGraph( H, G );
    boost::topological_sort(G, std::back_inserter(torder));

    #ifdef DEBUG
    undirGraphT Gp;
    boost::copy_graph(G,Gp);
    std::vector<int> component(num_vertices(Gp));
    int num = connected_components(Gp, &component[0]);
    LOG_INFO(log) << "** NUMBER OF CONNECTED COMPONENTS = " << num << "\n\n";
    #endif 

    // Remove all of the vertices from "order" if they are not in vset
    LOG_INFO(log) << "ORDER BEFORE = " << torder.size() << "\n";
    std::copy_if(torder.begin(), torder.end(), std::back_inserter(order), [ = ](const int & id) -> bool { return vset.find(id) != vset.end(); } );
    LOG_INFO(log) << "ORDER AFTER = " << order.size() << "\n";
    
    #ifdef DEBUG
    size_t start = 0; size_t end = H->order();
    for ( auto ind : boost::irange(start,end) ) {
        if ( vset.find(ind) == vset.end() ){
            auto k = H->vertex(ind);
            std::cerr  << "Lost {" << tree->getNodeName(k.u()) << ", " << tree->getNodeName(k.v()) <<
                          ", (" << k.f() << ", " << k.r() << ")} (" <<  
                          k.connectU() << ", " << k.connectV() << ") in the topological sort\n";
            if ( differentExtantNetworks(ti, k.u(), k.v()) ){
                std::cerr << "WE SHOULD KNOW BETTER; THE VERTICES ARE IN DIFFERENT EXTANT NETWORKS\n";
            }

            if ( ti.inSubnodesOf(k.u(), k.v()) ) { std::cerr << "U in subnodes of V\n"; }
            if ( ti.inSubnodesOf(k.v(), k.u()) ) { std::cerr << "V in subnodes of U\n"; }
            if ( isLost(k.u()) ) { std::cerr << "U lost\n"; }
            if ( isLost(k.v()) ) {std::cerr << "V lost\n"; }

        }
    }
    #endif
}

/**
 *  Compute the penalty for this edge to exist based on difference
 *  between the existence intervals of the endpoints and the
 *  penalty factor.
 */
template <typename T>
double existencePenalty( const TreeInfo &ti, const T &vert, double penalty, double travWeight ) {
    // Original (March 5)    
    // auto dist = ti.intervalDistance( vert.u(), vert.v() );
    // auto mult = (vert.state() == FlipState::none) ? 1.0 : 0.5;    
    // auto penatly = (travWeight > 0.0 and dist > 0.0) ? 1.0 + mult*std::exp(dist*penalty) : 1.0;
    // double res1 = travWeight * (( penalty > 1.0 ) ? penalty : 1.0);
    // return res1;

    // double height = 2.0 * 3.23;
    // auto d1 = const_cast<TreeInfo&>(ti).extantInterval[vert.u()].birth;
    // auto d2 = const_cast<TreeInfo&>(ti).extantInterval[vert.v()].birth;

    double res2 = 0.0;
    auto dist = ti.intervalDistance( vert.u(), vert.v() );
    auto mult = (vert.state() == FlipState::none) ? 1.0 : 1.0;    
    if ( dist > 0.0 ) { dist = 1.0; }
    //if ( dist > 0.5 ) { dist = 10.0; } else { dist = 0.0; }
    //res2 = travWeight  + mult * (std::exp(penalty*dist*dist) - 1.0);
    res2 = travWeight + travWeight * penalty;// * std::exp(dist);

    // if ( travWeight != res2 ) {
    //     std::cerr << "before penalty = " << travWeight << ", after penalty = " << res2 << "\n";
    // }

    return res2;

}


unique_ptr<ForwardHypergraph>  buildMLSolutionSpaceGraph( const TreePtrT &t,
        const TreeInfo &ti,
        Model &model,
        bool directed ) {

    cpplog::FileLogger log( "log.txt", true );

    boost::timer::auto_cpu_timer timer;

    auto isLeaf = [&]( const int & nid ) -> bool { return t->isLeaf(nid); };
    auto isInternal = [&]( const int & nid ) -> bool { return !t->isLeaf(nid); };
    //auto isLost = [&]( int nid ) -> bool { return (t->getNodeName(nid)).find("LOST") != std::string::npos; };
    auto isLost = [&]( int nid ) -> bool { return false; };


    auto vstr = [&]( const FlipKey & vert ) -> string {
        auto uname = t->getNodeName(vert.u());
        auto vname = t->getNodeName(vert.v());
        if (uname > vname) {
            auto tmp = vname;
            vname = uname;
            uname = tmp;
        }

        string fstr = "";
        if ( vert.state() == FlipState::both ) {
            fstr = "<-->";
        } else if ( vert.state() == FlipState::none ) {
            fstr = "X";
        } else {
            std::abort();
        }
        return "[" + uname + ", " + vname + "] : " + fstr;
    };


    /**
    * Determines if the transition from the key parent to the key child is valid.
    * Any transition in which the state of the interaction is not altered (not flipped) between the 
    * parent and child is valid.  If the state of the interaction is flipped, the transition is only
    * valid if the interacting nodes ( (u,v) in the parent ) exist in the same species.
    */
    auto isValidTransition = [&]( const FlipKey& parent, const FlipKey& child ) -> bool {
        return true;
        if ( parent.state() == FlipState::none and (parent.state() != child.state()) ) {
            auto us = dynamic_cast<bpp::BppString*>( t->getNodeProperty(parent.u(),"S") )->toSTL();
            auto vs = dynamic_cast<bpp::BppString*>( t->getNodeProperty(parent.v(),"S") )->toSTL();
            return us == vs;
        } 
        return true;        
    };


    auto rootId = t->getRootId();
    auto rootName = t->getNodeName( rootId );

    auto fauxRoot = "#preroot#";

    unique_ptr<ForwardHypergraph> slnSpaceGraph( new ForwardHypergraph() );

    auto addIncomingHyperedge = [ = , &slnSpaceGraph, &t] (const FlipKey & k ) {
        typedef std::vector<FlipKey> FlipKeyVec;

        auto addCartesianProduct = [&] ( 
            std::vector<FlipKeyVec> &iterKeys, 
            std::vector<FlipKey> &fixedKeys,
            std::function<bool(std::vector<FlipKey>&)>& acceptEdge ) -> void {

            std::vector<size_t> indexes( iterKeys.size(), 0 );
            bool done = false;
            while ( !done ) {
                FlipKeyVec res;
                res.reserve( iterKeys.size() );
                bool isValidEdge = true;
                for ( size_t i = 0; i < indexes.size(); ++i ) {
                    auto& childKey = iterKeys[i][indexes[i]];
                    if ( !isValidTransition(k, childKey) ) { 
                        isValidEdge = false; break; 
                    }
                    res.push_back(childKey);
                }
                res.insert( res.end(), fixedKeys.begin(), fixedKeys.end() );
                double w = 1.0;
                if ( isValidEdge and acceptEdge(res) ) { 
                    // Print the hyperedge
                    // std::cerr << "Added hyperedge " << vstr(k) << " => ";
                    // for ( auto tk : res ){ std::cerr << vstr(tk) << ", ";}
                    // std::cerr << "\n";
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

        auto dirKey = k.state();
        double canonicalDerivCost = 0.0;
        const int LN = k.u();
        const int RN = k.v();
        
        const int U = k.u();
        const int V = k.v();
        int LU, RU, LV, RV;
        LU = RU = LV = RV = -1;
        if (isInternal(U)) { 
            auto sons = t->getSonsId(U);
            LU = sons[0]; RU = sons[1];
            if ( LU > RU ) { std::swap(LU,RU); }
        }
        if (isInternal(V)) {
            auto sons = t->getSonsId(V);
            LV = sons[0]; RV = sons[1];
            if ( LV > RV ) { std::swap(LV,RV); }
        }

        typedef std::vector< FlipKey> EdgeTail;
        typedef std::function<bool(EdgeTail&)> CheckEdgeFunction;
            
        // Don't filter any extra edges
        CheckEdgeFunction edgeAlwaysOk = [&] ( std::vector< FlipKey >& edge ) -> bool { 
            return edge.size() > 0;
        };

        // Require at least one edge --- fw, rev, or both --- to exist
        CheckEdgeFunction atLeastOneEdge = [&] ( EdgeTail& edge ) -> bool {
            if ( edge.size() > 0 ) {
                for ( auto& k : edge ) {
                    if ( k.state() != FlipState::none ) { return true; }
                }
            }
            return false;
        };

        // Self loop
        if ( k.isSelfLoop() ) {
            if (isInternal(U)) {
                /** Possible states
                * {LU,LU}  {RU, RU}  {LU, RU}
                *    0        0         0
                *    0        0         1
                *    0        1         0
                *    0        1         1
                *    1        0         0
                *    1        0         1
                *    1        1         0
                *    1        1         1
                */
                std::vector< FlipKeyVec > iterKeys;
                std::vector<FlipKey> none;
                if ( !isLost(LU) ) {
                    iterKeys.push_back( {{LU, LU, FlipState::both},{LU, LU, FlipState::none}} );
                }
                if ( !isLost(RU) ) {
                    iterKeys.push_back( {{RU, RU, FlipState::both},{RU, RU, FlipState::none}} );
                }
                if ( not differentExtantNetworks(ti, LU, RU) and ! (isLost(LU) or isLost(RU)) ) {
                    iterKeys.push_back( {{LU, RU, FlipState::both},{LU, RU, FlipState::none}} );
                }
                addCartesianProduct(iterKeys, none, edgeAlwaysOk);
            
            } // isInternal(U)

        } else { // not a self-loop

            std::vector< FlipKeyVec > iterKeys;
            std::vector< FlipKey > none;
            /**  Covers the cases where U duplicates first
            *    s(LU,V)   s(RU,V)
            *    -------   -------
            *       T          T
            *       T          F
            *       F          T
            *       F          F
            */
            if ( isInternal(U) ) {
                if ( LU > RU ) { std::swap(LU, RU); }

                if ( not differentExtantNetworks(ti, LU, V) && !(isLost(LU) || isLost(V)) ) {
                    iterKeys.push_back( {FlipKey(LU, V, FlipState::both), 
                                           FlipKey(LU, V, FlipState::none)} ); 
                }

                if ( not differentExtantNetworks(ti, RU, V) && !(isLost(RU) || isLost(V)) ) {
                    iterKeys.push_back( {FlipKey(RU, V, FlipState::both), 
                                           FlipKey(RU, V, FlipState::none)} ); 
                }

            }
            addCartesianProduct(iterKeys, none, edgeAlwaysOk);
            iterKeys.clear();


            /**  Covers the cases where V duplicates first
            *    s(U,LV)   s(U,RV)
            *    -------   -------
            *       T          T
            *       T          F
            *       F          T
            *       F          F
            */
            if ( isInternal(V) ) {
                if ( LV > RV ) { std::swap(LV, RV); }

                if ( not differentExtantNetworks(ti, LV, U) && !(isLost(LV) || isLost(U)) ) {
                    iterKeys.push_back( {FlipKey(U, LV, FlipState::both),
                                           FlipKey(U, LV, FlipState::none)});
                }

                if ( not differentExtantNetworks(ti, RV, U) && !(isLost(RV) || isLost(U)) ) {
                    iterKeys.push_back( {FlipKey(U, RV, FlipState::both),
                                        FlipKey(U, RV, FlipState::none)});
                }

            }
            addCartesianProduct(iterKeys, none, edgeAlwaysOk);
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
    for ( uit = tbegin; uit != tend; ++uit ) {
        int u = *uit;
        for ( vit = uit; vit != tend; ++vit ) {
            int v = *vit;

            // Add the pair if they meet the following conditions
            if ( ( (u == v) && !isLost(u) ) || // they are the same node (self-loop)
                 (!differentExtantNetworks(ti, u, v) && // they have progeny in the same species
                 !(ti.inSubnodesOf(u, v) ||  ti.inSubnodesOf(v, u)) && // neither is a progeny of the other
                 !(isLost(u) || isLost(v)) )) { // neither of them is lost
                
                // Add a vertex representing no edge between them
                slnSpaceGraph->addVertex( FlipKey( u, v, FlipState::none ) );
                if ( ! t->isRoot(v) ) {
                    // Add a vertex representing an edge between them (or both edges in a dir. graph)
                    slnSpaceGraph->addVertex( FlipKey( u, v, FlipState::both ) );
                }
                if ( directed && u != v ) {
                    // If this is a directed reconstruction, add both the forward and reverse edges
                    slnSpaceGraph->addVertex( FlipKey( u, v, FlipState::forward ) );
                    slnSpaceGraph->addVertex( FlipKey( u, v, FlipState::reverse ) );
                }

            } // end node addition
        } // end vit
    } // end uit

    // The hypergraph has N vertices
    auto N = slnSpaceGraph->order();
    ez::ezETAProgressBar showProgress(N); showProgress.start();


    // For each vertex
    for ( size_t i = 0; i < N; ++i, ++showProgress ) {
        auto k = slnSpaceGraph->vertex(i);
        addIncomingHyperedge( k );
    }

    //cerr << "\n";
    LOG_INFO(log) << "Hypergraph size = " << slnSpaceGraph->size() << "\n";
    return slnSpaceGraph;
}

/**
 *  Add the appropriate hypervertices representing states
 *  between u and v.  This function generates hypervertices 
 *  that carry the extra information about the non-ambiguous 
 *  (i.e. only interaction state changes matter, and not the order of 
 *   dupication events leading to them) derivation.
 */
void addHyperNodeUnambiguous (
    int u, 
    int v, 
    FlipState fs,
    const TreePtrT& t,
    const TreeInfo &ti,
    unique_ptr<ForwardHypergraph>& slnSpaceGraph ) {

    if ( u > v ) { std::swap(u,v); }
    
    slnSpaceGraph->addVertex( FlipKey( u, v, fs, MustConnect::none) );
    if ( u != v ) { // connect states only matter for heterodimers
        if ( Utils::Trees::sameSpecies(t, u, v) ) {
            slnSpaceGraph->addVertex( FlipKey( u, v, fs, MustConnect::both) );
        }

        if ( !t->isLeaf(u) and Utils::Trees::isDescendantSpecies(ti, v, u) )  {// only add this node if u is internal
            slnSpaceGraph->addVertex( FlipKey( u, v, fs, MustConnect::right) );
        }

        if ( !t->isLeaf(v) and Utils::Trees::isDescendantSpecies(ti, u, v) )  { // only add this node if v is internal
            slnSpaceGraph->addVertex( FlipKey( u, v, fs, MustConnect::left) );
        }
    }

}

/**
 *  Add the appropriate hypervertices representing states
 *  between u and v.  This function generates hypervertices 
 *  that don't encode any information about the interaction 
 *  changes beneath them.  The hypergraph built with these
 *  vertices will count two different ordering of duplication
 *  events (even if they imply the same eventual set of interactions)
 *  as two distinct derivations.
 */
void addHyperNodeAmbiguous (
    int u, 
    int v, 
    FlipState fs,
    const TreePtrT& t,
    const TreeInfo &ti,
    unique_ptr<ForwardHypergraph>& slnSpaceGraph ) { 

    if ( u > v ) { std::swap(u,v); }
    slnSpaceGraph->addVertex( FlipKey( u, v, fs, MustConnect::none) );    
}

/**
* This function builds the hypergraph structure from the dynamic programming recurrences.
**/
unique_ptr<ForwardHypergraph>  buildSolutionSpaceGraph( 
    const TreePtrT &t,
    const TreeInfo &ti,
    double cc,
    double dc,
    double penalty,
    bool directed,
    MultiOpt::DerivationType dtype ) {

    typedef std::function<bool(const FlipKey&)> AddHyperEdgeFuncT;

    cpplog::FileLogger log( "log.txt", true );

    boost::timer::auto_cpu_timer timer("Building Hypergraph [%ws wall, %us user + %ss system = %ts CPU (%p%)]\n");

    std::unordered_set<int> leafIDs;
    std::unordered_set<int> lostIDs;

    for ( auto lid : t->getLeavesId() ) {
        if ( t->getNodeName(lid).find("LOST") != std::string::npos ) {
            lostIDs.insert(lid);
        }
        leafIDs.insert(lid);
    }

    auto isLeaf = [&]( const int & nid ) -> bool { return leafIDs.find(nid) != leafIDs.end(); };
    auto isInternal = [&]( const int & nid ) -> bool { return !isLeaf(nid); };
    auto isLost = [&]( int nid ) -> bool { return lostIDs.find(nid) != lostIDs.end(); }; //return (t->getNodeName(nid)).find("LOST") != std::string::npos; };

    auto rootId = t->getRootId();
    auto rootName = t->getNodeName( rootId );

    auto fauxRoot = "#preroot#";

    // Create the hypergraph we'll be dealing with
    unique_ptr<ForwardHypergraph> slnSpaceGraph( new ForwardHypergraph() );

    // The cost of transitions between different states and self-loops
    costMapT costMap( getCostDict(cc, dc, directed) );
    selfCostMapT selfLoopCostMap( getSelfLoopCostDict(cc, dc, directed) );

    /**
    * A term (hypervertex) is valid if neither u nor v is lost and the extant
    * networks descending from u and v overlap.
    */
    auto isValidTerm = [&]( const FlipKey& fk ) -> bool {
        return (!differentExtantNetworks(ti, fk.u(), fk.v()) && !(isLost(fk.u()) || isLost(fk.v())));
    };

    /**
    * Convenience function that performs the same test as 'isValidTerm' but which takes
    * a vector of node ids (assumed of length >=2) instead of a FlipKey.
    */
    auto isInvalidPair = [&]( const std::vector<int>& uv ) -> bool {
        return (differentExtantNetworks(ti, uv[0], uv[1]) || isLost(uv[0]) || isLost(uv[1]));
    };

    /**
    * Determines if the transition from the key parent to a child with the FlipState 'cs' is valid.
    * Any transition in which the state of the interaction is not altered (not flipped) between the 
    * parent and child is valid.  If the state of the interaction is flipped, the transition is only
    * valid if the interacting nodes ( (u,v) in the parent ) exist in the same species.
    */
    auto flipIsValid = [&]( const FlipKey& parent, FlipState cs ) {
        auto ps = parent.state();
        if ( ps != cs ) {
            auto us = dynamic_cast<bpp::BppString*>( t->getNodeProperty(parent.u(),"S") )->toSTL();
            auto vs = dynamic_cast<bpp::BppString*>( t->getNodeProperty(parent.v(),"S") )->toSTL();
            return us == vs;
        } 
        return true;                
    };

    /**
    * Convenience function that performs the same test as 'flipIsValid', but takes the 
    * parent and child keys directly as arguments.
    */
    auto isValidTransition = [&]( const FlipKey& parent, const FlipKey& child ) -> bool {
        return flipIsValid( parent, child.state() );
    };

    /**
    *  Print the key 'k' in a human readable format.
    */
    auto printKey = [&]( const FlipKey& k ) -> void {
        std::cerr << "{ " << t->getNodeName(k.u()) << ", " << t->getNodeName(k.v()) << ", (" <<
            k.f() << ", " << k.r() << ")} ( " << k.connectU() << ", " << k.connectV() << ")\n";
    };

    /**
    * Determines if 'connectState' constitutes a valid connection state between 'u' and 'v'.
    * The connection state of 'none' is always valid, 'left' is valid only if
    * 'v' is not a leaf and 'u' is a descendant species of 'v' (or the same species).
    * The 'right' connect state is analogous to 'left' but with 'u' and 'v' swapped, and
    * the 'both' connect state is only valid of 'u' and 'v' exist in the same species.
    */
    auto validConnectState = [&]( int u, int v, MustConnect connectState ) -> bool {
        switch ( connectState ) {
            case MustConnect::none:
              return true;
            case MustConnect::left:
              return (!t->isLeaf(v)) and Utils::Trees::isDescendantSpecies(ti, u, v);
            case MustConnect::right:
              return (!t->isLeaf(u)) and Utils::Trees::isDescendantSpecies(ti, v, u);
            case MustConnect::both:
              return Utils::Trees::sameSpecies(t,u,v);
        }
    };

    /**
    * Adds the appropriate incoming hyperedges to the hypervertex 'k'.
    */
    AddHyperEdgeFuncT addIncomingHyperedgeUnambiguous = [ = , &costMap, &selfLoopCostMap, &slnSpaceGraph, &penalty, &t] (
        const FlipKey & k //!< The vertex to which we'll add the incoming edges
        ) -> bool {
        
        /**
        * Adds all possible (and valid) combinations of edges defined by 
        * the hypervertex set 'targetNodes' and the connectOptions set 'connectOptions'
        * to the with a head hypervertex of 'k'.
        */
        auto addPossibleEdges = [&]( 
            std::vector< std::vector<int> >& targetNodes, 
            std::vector< MustConnect >& connectOptions, 
            const FlipKey& k,
            FlipState fs,
            double cost, 
            std::function<bool(std::vector<FlipKey>&)>& acceptEdge ) -> bool {

            using boost::next_partial_permutation;
            
            bool added{false};

            try {
                // If we're trying to perform a flip (the flip state of the parent and potential
                // children is different) but the parent nodes are in different species, then none
                // of the potential edges is valid
                if ( !flipIsValid(k, fs) ) { return false; }

                // Remove any invalid vertices in the tail (i.e. vertices whose constituent nodes
                // are lost or are in different extant networks.).
                auto newEnd = std::remove_if(targetNodes.begin(), targetNodes.end(), isInvalidPair);
                targetNodes.erase(newEnd, targetNodes.end());
                if (targetNodes.size() == 0) { return false; } //std::cerr << "invalid edge\n"; return false;}
                
                int numConnectOptions = connectOptions.size();
                int numTargets = targetNodes.size();

                std::vector<size_t> inds;
                inds.reserve( numConnectOptions * numTargets );

                for( auto i : boost::irange(0, numConnectOptions) ) {
                    for( auto j : boost::irange(0, numTargets) ) { inds.push_back(i); }
                }

                std::vector< FlipKey > edge;
                edge.reserve(numTargets);

                std::unordered_set< std::string > uniqueEdges;

                size_t numAdded = 0;
                do {

                  // Build up each possible edge by pairing all of the target nodes with
                  // each permutation of connect states
                  for( auto tn : boost::irange(0,numTargets) ) {  // Build "edge"

                    // The pair (u,v) is a current tail node of "edge"
                    auto u = targetNodes[tn][0]; auto v = targetNodes[tn][1];                    
                    // The connect state we'll assign to (u,v)                    
                    auto mustConnect = connectOptions[inds[tn]];
                    assert(u<=v);
                    if ( u > v ) { std::cerr << "swapping\n"; std::swap(u,v); swapEndpoints(mustConnect); }

                    // The actual key representing this tail node of "edge".
                    FlipKey fk( u, v, fs, mustConnect);
                    
                    // If the connect state is valid, then add this node to the edge, otherwise
                    // this edge does not represent a valid recurrence, so skip it.
                    if ( validConnectState(u,v,mustConnect) ) { 
                        edge.push_back(fk);
                    } else {
                        edge.clear(); break;
                    }

                  }

                  if ( acceptEdge(edge) ) {
                      added = true;
                      numAdded += 1;
                      slnSpaceGraph->addEdge( edge, k, cost );
                  }              
                  edge.clear();

                } while ( next_partial_permutation( inds.begin(), inds.begin()+numTargets, inds.end() ) );
                  
           } catch ( FlipKey& fk ) {

                std::cerr << "Error adding edge; couldn't find node : ";
                std::cerr << "{ " << t->getNodeName(fk.u()) << ", " << t->getNodeName(fk.v()) << ", (" <<
                    fk.f() << ", " << fk.r() << ") } (" << 
                    fk.connectU() << ", " << fk.connectV() << "), ";

           } catch ( ... ) {
             std::cerr << "failed to add an edge!\n";
           }
           return added;
        };    

        bool addedEdge{false};

        auto f = k.f(); auto r = k.r();        
        auto state = k.state();
        auto flippedState = fsBoth[state];

        MustConnect mustConnect{ k.mustConnect() };

        bool connectU{k.connectU()};
        bool connectV{k.connectV()};

        const int U = k.u();
        const int V = k.v();
        auto UIsInternal = isInternal(U);
        auto VIsInternal = isInternal(V);

        int LU, RU, LV, RV;
        LU = RU = LV = RV = -1;

        if ( UIsInternal ) { 
            auto sons = t->getSonsId(U);
            assert(sons.size() == 2);
            LU = sons[0]; RU = sons[1]; 
            if ( LU > RU ) { std::swap(LU,RU); }
        }
        if ( VIsInternal ) { 
            auto sons = t->getSonsId(V);
            assert(sons.size() == 2);
            LV = sons[0]; RV = sons[1]; 
            if ( LV > RV ) { std::swap(LV,RV); }
        }

        // Handle the case where this vertex denotes a self loop
        if (k.arity() == 1) {

            // There are only incoming edges to add if this vertex isn't a leaf
            if (isInternal(U)) {

                std::vector< MustConnect > lrStates = { 
                    MustConnect::none, MustConnect::both, MustConnect::right, MustConnect::left
                };

                // 1 -- we don't flip the self-loop
                assert( mustConnect == MustConnect::none );
                auto noFlipLU = FlipKey( LU, LU, state, mustConnect );
                auto noFlipRU = FlipKey( RU, RU, state, mustConnect );
                // 2 -- we flip the self loop
                auto dualFlipLU = flipBoth( noFlipLU );
                auto dualFlipRU = flipBoth( noFlipRU );

                vector<FlipKey> noFlipEdge = vector<FlipKey>();
                vector<FlipKey> dualFlipEdge = vector<FlipKey>();

                bool empty{true};
                // we can only recurse into U if it's not lost
                if ( !isLost(LU) ) {
                    noFlipEdge.push_back( noFlipLU );
                    dualFlipEdge.push_back( dualFlipLU );
                    empty = false;
                }
                // same with V
                if ( !isLost(RU) ) {
                    noFlipEdge.push_back( noFlipRU );
                    dualFlipEdge.push_back( dualFlipRU );
                    empty = false;
                }

                // We only consider the state between U & V if they are part of the same
                // extant network and neither is lost
                if ( !differentExtantNetworks(ti, LU, RU) and !(isLost(LU) or isLost(RU)) ) {

                    // Pair up the homodimer vertices with each of the possible connect states 
                    // of the heterodimer vertices
                    for ( auto& cstate : lrStates ) {

                        if ( validConnectState(LU, RU, cstate) ) {
                            noFlipEdge.push_back( FlipKey(LU, RU, state, cstate ));
                            slnSpaceGraph->addEdge( noFlipEdge, k, 0.0 );

                            auto ckey = FlipKey(LU, RU, flippedState, cstate );
                            
                            // if ( isValidTransition(k, ckey) ) { 
                            // DON'T NEED THIS TEST HERE; The same node should ALWAYS be from the same species
                            
                            dualFlipEdge.push_back( ckey );
                            auto w = round3(get<0>(selfLoopCostMap[ k.state() ][ (!f || !r) ]));
                            w = existencePenalty(ti, k, penalty, w);

                            if ( std::isfinite(w) ) {
                                addedEdge = true;
                                slnSpaceGraph->addEdge( dualFlipEdge, k, w );
                            }

                            // get rid of the heterodimer so we can make the next edge
                            noFlipEdge.pop_back();
                            dualFlipEdge.pop_back();
                        }
                    }

                } else if ( !empty ) { // If we don't need to consider edges between LU & RU
                    slnSpaceGraph->addEdge( noFlipEdge, k, 0.0 );

                    auto w = round3(get<0>(selfLoopCostMap[ k.state() ][ (!f || !r) ]));
                    w = existencePenalty(ti, k, penalty, w);
                    if ( std::isfinite(w) ) {
                        addedEdge = true;
                        slnSpaceGraph->addEdge( dualFlipEdge, k, w );
                    }
                }

            }

        } else { // if U and V are different --- this isn't a self-loop

          typedef std::vector<FlipKey> EdgeTail;
          typedef std::function<bool(EdgeTail&)> CheckEdgeFunction;
            
          // Don't filter any extra edges
          CheckEdgeFunction edgeAlwaysOk = [&] ( std::vector< FlipKey >& edge ) -> bool { 
            return edge.size() > 0;
          };
            
          // For at least one flip key in this edge, the state of connectU is true
          CheckEdgeFunction childrenMustConnectU = [&] ( std::vector< FlipKey >& edge ) -> bool {
            if ( edge.size() > 0 ) {
                for ( auto& k : edge ){
                    if ( k.connectU() ) { return true; }
                }
            }
            return false;
          };

        // For at least one flip key in this edge, the state of connectV is true
        CheckEdgeFunction childrenMustConnectV = [&] ( std::vector< FlipKey >& edge ) -> bool { 
            if ( edge.size() > 0 ) {
                for ( auto& k : edge ){
                    if ( k.connectV() ) { return true; }
                }
            }
            return false;
        };            

        if ( mustConnect == MustConnect::both ) {
          if ( UIsInternal or VIsInternal ) {
            
            // Since connectU  and connectV are both true, then something
            // must flip the state to U and something must flip the state to V.
            // The only way this is possible without creating a blocking loop is by 
            // flipping {U,V}.  Thus, we consider all possible ways of recursing
            // on U and V, but we only consider recursions where we flip the state
            // of the interaction.
            
            // The nodes into which we'll recurse
            std::vector< std::vector<int> > targetNodes;
            // The set of connect options we'll allow in these nodes
            std::vector< MustConnect > connectOptions{ 
                MustConnect::none, MustConnect::left, MustConnect::right, MustConnect::both 
            }; 


            // The cost of performing the flip
            auto w = round3(get<0>(costMap[ state ][ flippedState ]));
            w = existencePenalty(ti, k, penalty, w);

            // If both U and V are internal we consider the cases of recursing into both
            // of them simultaneously.  Any path which takes one of these edges will consider
            // no further flips to U or V.
            if ( UIsInternal and VIsInternal ) { 
                targetNodes = { {LU,LV}, {LU,RV}, {RU,LV}, {RU,RV} };
                assert( LU <= LV ); assert( LU <= RV );
                assert( RU <= LV ); assert( RU <= RV );
                addedEdge |= addPossibleEdges( targetNodes, connectOptions, k,  flippedState, w, edgeAlwaysOk );
            }

            // If V is internal, then we can recurse into its children.  However,
            // to avoid ambiguity, we will only consider the cases where a descendant of V
            // still flips the state to U (otherwise, we only consider recursing on both
            // which is handled above).
            if ( VIsInternal ) {
               targetNodes = { {U,LV}, {U,RV} };
               assert(U <= LV); assert(U <= RV);
               addedEdge |= addPossibleEdges( targetNodes, connectOptions, k, flippedState, w, childrenMustConnectU );
            }

            // If U is internal, then we can recurse into its children.  However,
            // to avoid ambiguity, we will only consider the cases where a descendant of U
            // still flips the state to V (otherwise, we only consider recursing on both
            // which is handled above).
            if ( UIsInternal ) {
                targetNodes = { {LU,V}, {RU,V} };
                assert(LU <= V); assert(RU <= V);
                addedEdge |= addPossibleEdges( targetNodes, connectOptions, k, flippedState, w, childrenMustConnectV );
            }

        } 
    } else if ( mustConnect == MustConnect::left ) {
        
        // Since connectV is 0, we *can not* flip {U,V}, but _something_ beneath V *must* still 
        // flip the state to U.
        // This implies that we only consider recursing on V and recursing without flipping.
        // Also, since something beneath V must still flip the state to U, we only consider
        // recursing to nodes with the flip state {1,0} or {1,1}
        
        if ( VIsInternal ) { // If we can't recurse on V then we're hosed
           assert(U < LV);
           assert(U < RV);
           std::vector< std::vector<int> > targetNodes{ {U,LV}, {U,RV} };
           std::vector< MustConnect > connectOptions{ 
            MustConnect::both, MustConnect::left, MustConnect::right, MustConnect::none 
           };
           addedEdge = addPossibleEdges( targetNodes, connectOptions, k, state, 0.0, childrenMustConnectU );
        } else { 
            std::cerr << "ERROR\n";
            std::abort();
        }                

    } else if ( mustConnect == MustConnect::right ) {
        
        // Since connectU is 0, we *can not* flip {U,V}, but _something_ beneath U *must* still 
        // flip the state to V.
        // This implies that we only consider recursing on U and recursing without flipping.
        // Also, since something beneath U must still flip the state to V, we only consider
        // recursing to nodes with the flip state {0,1} or {1,1}

        if ( UIsInternal ) { // If we can't recurse on U then we're hosed
            assert(LU < V);
            assert(RU < V);
            std::vector< std::vector<int> > targetNodes = { {LU,V}, {RU,V} };
            std::vector< MustConnect > connectOptions{ 
                MustConnect::both, MustConnect::left, MustConnect::right, MustConnect::none 
            };
            
            addedEdge = addPossibleEdges( targetNodes, connectOptions, k, state, 0.0, childrenMustConnectV );
        } else { 
            std::cerr << "ERROR\n";
            std::abort();
        }

    } else if ( mustConnect == MustConnect::none ) {
        // If neither node is internal, then we can't recurse at all
        if ( UIsInternal or VIsInternal )  {

            // We're not allowed to change the state to U or V, so we must not flip.

            // The nodes into which we'll recurse
            std::vector< std::vector<int> > targetNodes;
            // The set of connect options we'll allow in these nodes
            std::vector< MustConnect > connectOptions;
            
            // If both U and V are internal
            if ( UIsInternal and VIsInternal ) { 
                targetNodes = { {LU,LV}, {LU,RV}, {RU,LV}, {RU,RV} };
                connectOptions = { MustConnect::both, MustConnect::left, MustConnect::right, MustConnect::none };
            } else

            // If U is a leaf, then the connect state of the childern must be {0,0}
            // since there is nothing below U to connect to V and since the connect
            // state at U,V is {0,0}, nothing below V is allowed to connect to U.
            if ( (!UIsInternal) and VIsInternal) {
                targetNodes = { {U,LV}, {U,RV} };
                connectOptions = { MustConnect::none };
            } else

            // If V is a leaf, the same argument as above applies and the connect stae
            // of all children must be {0,0}.
            if ( UIsInternal and (!VIsInternal) ) {
                targetNodes = { {LU,V}, {RU,V} };
                connectOptions = { MustConnect::none };                  
            }
            addedEdge = addPossibleEdges( targetNodes, connectOptions, k, state, 0.0, edgeAlwaysOk );            
        }

    }
        }
        return addedEdge;
    };

    /**
    * Adds the appropriate incoming hyperedges to the hypervertex 'k'.
    */
    AddHyperEdgeFuncT addIncomingHyperedgeAmbiguous = [ = , &costMap, &selfLoopCostMap, &slnSpaceGraph, &penalty, &t] (
        const FlipKey & k //!< The vertex to which we'll add the incoming edges
        ) -> bool {


            auto state = k.state();
            auto flippedState = fsBoth[state];
            auto f = k.f(); auto r = k.r();

            double canonicalDerivCost = 0.0;
            const int U = k.u();
            const int V = k.v();
            auto UIsInternal = isInternal(U);
            auto VIsInternal = isInternal(V);

            int LU, RU, LV, RV;
            LU = RU = LV = RV = -1;

            if ( UIsInternal ) { 
                auto sons = t->getSonsId(U);
                LU = sons[0]; RU = sons[1]; 
                if ( LU > RU ) { std::swap(LU,RU); }
            }
            if ( VIsInternal ) { 
                auto sons = t->getSonsId(V);
                LV = sons[0]; RV = sons[1]; 
                if ( LV > RV ) { std::swap(LV,RV); }
            }

            if ( k.arity() == 1 ) { // self-loop

                // There are only incoming edges to add if this vertex isn't a leaf
                if (isInternal(U)) {
                  // 1 -- we don't flip the self-loop
                  auto noFlipLU = FlipKey( LU, LU, state );
                  auto noFlipRU = FlipKey( RU, RU, state );
                  // 2 -- we flip the self loop
                  auto dualFlipLU = flipBoth( noFlipLU );
                  auto dualFlipRU = flipBoth( noFlipRU );
                  
                  vector<FlipKey> noFlipEdge;
                  vector<FlipKey> dualFlipEdge;
                  
                  bool empty{true};
                  // we can only recurse into LU if it's not lost
                  if ( !isLost(LU) ) {
                    noFlipEdge.push_back( noFlipLU );
                    dualFlipEdge.push_back( dualFlipLU );
                    empty = false;
                  }
                  // same with RU
                  if ( !isLost(RU) ) {
                    noFlipEdge.push_back( noFlipRU );
                    dualFlipEdge.push_back( dualFlipRU );
                    empty = false;
                  }

                  // We only consider the state between U & V if they are part of the same
                  // extant network and neither is lost
                  if ( !differentExtantNetworks(ti, LU, RU) and !(isLost(LU) or isLost(RU)) ) {
                    noFlipEdge.push_back( FlipKey(LU, RU, state ));
                    auto flipBetween = FlipKey(LU, RU, flippedState );
                    if ( isValidTransition(k, flipBetween) ) {
                        dualFlipEdge.push_back( flipBetween );
                    }
                    empty = false; 
                  }

                  if ( !empty ) { // if there's something to add

                    slnSpaceGraph->addEdge( noFlipEdge, k, 0.0 );
                    auto opState = ( k.state() == FlipState:: none ) ? true : false;
                    auto w = round3(get<0>(selfLoopCostMap[ k.state() ][ opState ]));
                    w = existencePenalty(ti, k, penalty, w);

                    if ( std::isfinite(w) ) {
                        slnSpaceGraph->addEdge( dualFlipEdge, k, w );
                    }
                 }

             } // isInternal(U)

            } else { // non-self-loop

              if ( isInternal(U) ) { // Recurse into U

                 // Don't flip
                 vector<FlipKey> noFlip = vector<FlipKey>();
                 auto noFlipLU = ( V <= LU ) ? FlipKey( V, LU, state) : FlipKey( LU, V, state);
                 auto noFlipRU = ( V <= RU ) ? FlipKey( V, RU, state) : FlipKey( RU, V, state);
                 if ( isValidTerm(noFlipLU) ) { noFlip.push_back(noFlipLU); }
                 if ( isValidTerm(noFlipRU) ) { noFlip.push_back(noFlipRU); }

                 // Flip
                 vector<FlipKey> dualFlip = vector<FlipKey>();
                 auto dualFlipLU = flipBoth(noFlipLU);
                 auto dualFlipRU = flipBoth(noFlipRU);
                 if ( isValidTerm(dualFlipLU) and isValidTransition(k, dualFlipLU)) { dualFlip.push_back(dualFlipLU); }
                 if ( isValidTerm(dualFlipRU) and isValidTransition(k, dualFlipRU)) { dualFlip.push_back(dualFlipRU); }

                if ( noFlip.size() > 0 ) {
                    slnSpaceGraph->addEdge( noFlip, k, canonicalDerivCost);
                }

                if ( dualFlip.size() > 0 ) {
                    auto w = round3(get<0>(costMap[ state ][ flippedState ])); 
                    w = existencePenalty(ti, k, penalty, w);

                    if ( std::isfinite(w) ) {
                        slnSpaceGraph->addEdge( dualFlip, k, w );
                    }
                }

              } // isInternal(U)

              if ( isInternal(V) ) { // Recurse into V

                 // Don't flip
                 vector<FlipKey> noFlip = vector<FlipKey>();
                 auto noFlipLV = ( U <= LV ) ? FlipKey( U, LV, state) : FlipKey( LV, U, state);
                 auto noFlipRV = ( U <= RV ) ? FlipKey( U, RV, state) : FlipKey( RV, U, state);
                 if ( isValidTerm(noFlipLV) ) { noFlip.push_back(noFlipLV); }
                 if ( isValidTerm(noFlipRV) ) { noFlip.push_back(noFlipRV); }

                 // Flip
                 vector<FlipKey> dualFlip = vector<FlipKey>();
                 auto dualFlipLV = flipBoth(noFlipLV);
                 auto dualFlipRV = flipBoth(noFlipRV);

                 if ( isValidTerm(dualFlipLV) and isValidTransition(k, dualFlipLV)) { dualFlip.push_back(dualFlipLV); }
                 if ( isValidTerm(dualFlipRV) and isValidTransition(k, dualFlipRV)) { dualFlip.push_back(dualFlipRV); }

                if ( noFlip.size() > 0 ) {
                    slnSpaceGraph->addEdge( noFlip, k, canonicalDerivCost);
                }
                if ( dualFlip.size() > 0 ) {
                    auto w = round3(get<0>(costMap[ state ][ flippedState ])); 
                    w = existencePenalty(ti, k, penalty, w);

                    if ( std::isfinite(w) ) {
                        slnSpaceGraph->addEdge( dualFlip, k, w );
                    }
                }

              } // isInternal(V)

          } // non-self-loop
    };

    // Set the function to add the hypervertices and hyperedges based on
    // what type of derivation we consider 'unique'.
    auto addHypernode = (dtype == MultiOpt::DerivationType::AllHistories) ? \
                         addHyperNodeAmbiguous : addHyperNodeUnambiguous;

    auto addIncomingHyperedge = (dtype == MultiOpt::DerivationType::AllHistories) ? \
                                 addIncomingHyperedgeAmbiguous : addIncomingHyperedgeUnambiguous;
    
    // Create the vector of all nodes
    auto nodes = t->getNodesId();

    auto tbegin = nodes.cbegin();
    auto tend = nodes.cend();
    auto uit = tbegin, vit = tbegin;

    // Loop through all pairs of proteins to determine
    // which pairs represent valid hypervertices in the
    // hypergraph.  Add the associated hypervertices for these
    // protein pairs.
    for ( uit = tbegin; uit != tend; uit++ ) {
        int u = *uit;
        
        for ( vit = uit; vit != tend; vit++ ) {
            int v = *vit;

            if ( ((u == v) and !isLost(u) ) or
                !(differentExtantNetworks(ti, u, v) or
                  ti.inSubnodesOf(u, v) or
                  ti.inSubnodesOf(v, u) or
                  isLost(u) or isLost(v)) ) {

                addHypernode(u, v, FlipState::none, t, ti, slnSpaceGraph);

                if ( ! t->isRoot(v) ) {
                    addHypernode(u, v, FlipState::both, t, ti, slnSpaceGraph);
                }
                if ( directed and u != v ) {
                    addHypernode(u, v, FlipState::forward, t, ti, slnSpaceGraph);
                    addHypernode(u, v, FlipState::reverse, t, ti, slnSpaceGraph);
                }
            }

        }
    }

    auto N = slnSpaceGraph->order();
    cerr << "Hypergraph order: " << slnSpaceGraph->order() << "\n";
    //ProgressDisplay showProgress(N);
    ez::ezETAProgressBar eta(N); eta.start();
    for ( size_t i = 0; i < N; ++i, ++eta ) {
        bool addedEdge{false};
        auto k = slnSpaceGraph->vertex(i);
        addedEdge = addIncomingHyperedge( k ); 
    }
    //eta.n = N;

    LOG_INFO(log) << "Hyergraph order = " << slnSpaceGraph->order() << "\n";
    LOG_INFO(log) << "Hypergraph size = " << slnSpaceGraph->size() << "\n";

    return slnSpaceGraph;
}

template< typename GT >
void MLLeafCostDict( 
    unique_ptr<ForwardHypergraph> &H, 
    TreePtrT &T, 
    GT &G, 
    bool directed, 
    double cc, 
    double dc, 
    slnDictT &slnDict 
    ) {
    /*
      Given the duplication tree T, the root vertex rv, the extant graph G and
      the constraints, fill in the values for the leaf nodes of the hypergraph
    */
    typedef typename boost::graph_traits< GT >::edge_descriptor EdgeT;
    boost::timer::auto_cpu_timer timer;
    auto undirected = !directed;
    auto isLost = [&]( int nid ) -> bool { return (T->getNodeName(nid)).find("LOST") != std::string::npos; };
    auto isLeaf = [&]( int nid ) -> bool { return  T->isLeaf(nid); };
    // Cost of going from the state encoded by a node to the state of a true graph edge

    auto none = FlipState::none;
    auto fw = FlipState::forward;
    auto rev = FlipState::reverse;
    auto both = FlipState::both;

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
    auto contains = [] ( const NodeSetT & s, NodeT e ) { return s.find(e) != s.end(); };

    NodeSetT extantNodes;

    Google<int>::Set leafIds;
    leafIds.set_empty_key(-1);

    for ( auto l : T->getLeavesId() ) { leafIds.insert(l); }

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

        assert( isLeaf(nd.u()) and isLeaf(nd.v()) );

        // Check to see if u, v, or both have been lost
        auto endOfMap = idToVertMap.end();
        bool lostU = ( idToVertMap.find( nd.u() ) == endOfMap ) || ( ! contains(extantNodes, idToVertMap[nd.u()]) );
        bool lostV = ( idToVertMap.find( nd.v() ) == endOfMap ) || ( ! contains(extantNodes, idToVertMap[nd.v()]) );

        /*
        if ( !( isLeaf(nd.u()) and isLeaf(nd.v()) ) ) {
            std::cerr << nameU << ", " << nameV << ", cost = 0.5\n";
            //std::abort();
            auto cost = 0.5;
            vector<size_t> ev;
            Google<Derivation::flipT>::Set es;
            es.set_empty_key( make_tuple(-1, -1, ""));
            nlost += 1;
            slnDict[n] = { {0, Derivation(cost, n, ev, es)} };
        }*/

        // The cost to / between lost nodes is always 0
        if ( lostU || lostV ) {
            //double lostCost = ( nd.state() == FlipState::none ) ? 0.91 : 0.09;
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


/**
*  Given the duplication tree T, the root vertex rv, the extant graph G and
*  the constraints, fill in the values for the leaf hypervertices of the hypergraph.
*  The cost of each leaf hypervertex can be determined based on the state of edges in the
*  extant network and the cost of the various transitions.
*/
template< typename GT >
void leafCostDict( 
    unique_ptr<ForwardHypergraph> &H, //!< Problem Hypergraph 
    TreePtrT &T,                      //!< Phylogeny
    TreeInfo& ti,                     //!< Extra tree information
    GT &G,                            //!< Extant graph
    bool directed,                    //!< Is the graph directed?
    double cc,                        //!< Edge creation cost
    double dc,                        //!< Edge deletion cost
    slnDictT &slnDict                 //!< Dictionary of subproblem solutions
    ) {
    
    typedef typename boost::graph_traits< GT >::edge_descriptor EdgeT;    
    typedef typename GT::vertex_descriptor NodeT;
    typedef unordered_set<NodeT> NodeSetT;

    boost::timer::auto_cpu_timer timer("Populating Base Cases [%ws wall, %us user + %ss system = %ts CPU (%p%)]\n");
    auto undirected = !directed;
    auto isLost = [&]( int nid ) -> bool { return (T->getNodeName(nid)).find("LOST") != std::string::npos; };

    auto printKey = [&]( const FlipKey& k ) -> void {
        std::cerr << "{ " << T->getNodeName(k.u()) << ", " << T->getNodeName(k.v()) << ", (" <<
            k.f() << ", " << k.r() << ")} ( " << k.connectU() << ", " << k.connectV() << ")\n";
    };

    auto none = FlipState::none;
    auto fw = FlipState::forward;
    auto rev = FlipState::reverse;
    auto both = FlipState::both;

    auto costDict = getCostDict(cc, dc, directed);
    auto selfLoopCostDict = getSelfLoopCostDict(cc, dc, directed);
    auto costFunDict = getCostFunDict(cc, dc, directed);
    auto selfLoopCostFunDict = getSelfLoopCostFunDict(cc, dc, directed );

    auto N = H->order();
    auto M = H->size();

    using boost::counting_range;

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

    // A hypernode is a leaf if it has no incoming edges
    for ( auto i : counting_range(size_t{0},N) ) {
        if (H->incident(i).size() == 0) { leafHypernodes.push_back(i); }
    }

    // Is the node e contained in the set s?
    auto contains = [] ( const NodeSetT & s, NodeT e ) -> bool {
        return s.find(e) != s.end();
    };

    NodeSetT extantNodes;
    unordered_map<int, NodeT> idToVertMap;
    
    // Iterate over all vertices in the extant graph
    auto vp = boost::vertices(G);
    for ( auto it = vp.first; it != vp.second; ++it ) {
        auto v = *it;
        auto idx = G[v].idx;
        // If the vertex is a leaf, than it's extant
        if ( T->isLeaf(idx) ) { extantNodes.insert(v); }
        idToVertMap[idx] = v;
    }

    size_t nlost = 0;
    size_t nef = 0;
    double tweight = 0.0;

    vector<size_t> ev;
    Google<Derivation::flipT>::Set effectiveEdges;
    effectiveEdges.set_empty_key( make_tuple(-1, -1, ""));

    auto setLeafWeight = [&]( costRepT& costFlipProb, FlipKey& nd, int n ) -> void {
        auto cost = round3(get<0>(costFlipProb));
        auto flip = get<1>(costFlipProb);
        if ( flip != "n" ) { effectiveEdges.insert( make_tuple(nd.u(), nd.v(), flip) ); }
        slnDict[n] = { {0, Derivation(cost, n, ev, effectiveEdges)} };
    };

    for ( auto n : leafHypernodes ) { // For every leaf hypernode
    
        auto nd = H->vertex(n);
        auto nameU = T->getNodeName(nd.u());
        auto nameV = T->getNodeName(nd.v());

        // Check to see if u, v, or both have been lost
        auto endOfMap = idToVertMap.end();
        bool lostU = isLost(nd.u());
        bool lostV = isLost(nd.v());

        // If either of the vertices contained in this hypernode
        // is not a leaf, then the
        if ( !(T->isLeaf(nd.u()) and T->isLeaf(nd.v())) ) {
            assert(lostU or lostV);
            vector<size_t> ev;
            Google<Derivation::flipT>::Set es;
            es.set_empty_key( make_tuple(-1, -1, ""));            
            //auto cost = std::numeric_limits<double>::infinity();
            auto cost = 0.0;
            slnDict[n] = { {0, Derivation( cost, n, ev, es)} };
            continue;
        }

        // The cost to / between lost nodes is always 0
        if ( lostU or lostV ) {
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
            auto f = nd.f();  // forward edge
            auto r = nd.r();  // reverse edge

            if ( not nd.isSelfLoop() ) { // If this isn't a self-loop
                EdgeT fedge, redge;
                bool d_f, d_r;
                tie(fedge, d_f) = edge(u, v, G);
                double w_f = d_f ? G[fedge].weight : 0.0;
                tie(redge, d_r) = edge(v, u, G);
                double w_r = d_r ? G[redge].weight : 0.0;
                #ifdef DEBUG
                if ( undirected ) {
                    assert( w_f == w_r );
                }
                #endif

                tweight += (w_f + w_r) * (d_f + d_r) * (nd.f() + nd.r());

                auto costFlipProb = costDict[ nd.state() ][ directionsToFlipState(d_f,d_r) ];
                //auto costFlipProb = costFunDict[ nd.state() ][ directionsToFlipState(d_f,d_r) ](w_f, w_r);
                setLeafWeight( costFlipProb, nd, n);
            
            } else { // if this is a self-loop
            
                assert( u == v );
                EdgeT e;
                bool hasSelfLoop;
                tie(e, hasSelfLoop) = edge(u, u, G);
                double w_l = hasSelfLoop ? G[e].weight : 0.0;

                tweight += w_l * (hasSelfLoop) * (nd.f() + nd.r());
                auto costFlipProb = selfLoopCostDict[ nd.state() ][ hasSelfLoop ];
                //auto costFlipProb = selfLoopCostFunDict[ nd.state() ][ hasSelfLoop ](w_l);                
                setLeafWeight( costFlipProb, nd, n);
            } // ( u != v )
        } // ( lostU || lostV )
    } // loop over leaf hypernodes
}

template <typename CostClassT>
tuple<double, BigInt> getCostCount( vector<vector<CostClassT>> &tkd,
                                  const vector<size_t, StackAllocator<size_t>> &bp,
                                  const size_t &eid,
                                  unique_ptr<ForwardHypergraph> &H ) {

    auto edge = H->edge(eid);
    double cost = edge.weight();
    BigInt count(1);
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
                const vector<size_t, StackAllocator<size_t>> &bp,
                const size_t &eid,
                unique_ptr<ForwardHypergraph> &H ) {
    auto edge = H->edge(eid);
    double cost = edge.weight();
    size_t i = 0;
    for ( auto tailNode : edge.tail() ) {
        // the cost class index
        size_t cci = bp[i];
        // We sum the costs and multiply the counts
        cost += tkd[tailNode][cci].cost();
        i++;
    }
    return round3(cost);
}

template <typename CostClassT>
BigInt getCount( vector< vector<CostClassT> > &tkd,
               const vector<size_t, StackAllocator<size_t>> &bp,
               const size_t &eid,
               unique_ptr<ForwardHypergraph> &H ) {
    auto edge = H->edge(eid);
    BigInt count(1);
    size_t i = 0;
    for ( auto tailNode : edge.tail() ) {
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
// (1, 3200)
//
vector< tuple<double, BigInt> > countEdgeSolutions(
    const double &ecost,               //!< Cost of the edge that we're counting
    const vector<size_t> &tailNodes,   //!< The set of tail nodes of this edge
    countDictT &countDict,             //!< The dictionary holding the # of slns for nodes
    const size_t &k,                   //!< The # of cost classes to consider
    bool printMe,                      //!< Print extra info
    unique_ptr<ForwardHypergraph> &H,  //!< The hypergraph
    TreePtrT &t                        //!< The phylogeny
    ) {

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
    std::vector< size_t, StackAllocator<size_t>> elemSizes;
    elemSizes.reserve(tailNodes.size());
    double cost = ecost;
    for ( const auto & t : tailNodes ) {
        elemSizes.push_back( countDict[t].size() );
        cost += countDict[t].front().score;
    }

    vector<dvsT> pq(1, make_tuple(cost, vector<size_t, StackAllocator<size_t>>(tailNodes.size(), 0)));
    QueueCmp<dvsT> ord;


    std::function< double( const vector<size_t, StackAllocator<size_t>>& ) > computeScore = [&] ( 
        const vector<size_t, StackAllocator<size_t>> &inds ) -> double {
        size_t numNodes = tailNodes.size();
        double cost = ecost;
        for ( size_t i = 0; i < numNodes; ++i ) {
            cost += countDict[ tailNodes[i] ][ inds[i] ].score;
        }
        return cost;
    };

    typedef tuple<double, BigInt> ccT;
    vector< ccT > edgeSlns;
    double epsilon = 1.00;//5e-1;
    /*
    double minElement = round3(cost - 0.5*epsilon);
    vector< ccT > edgeSlns;
    for( size_t i = 1; i < k+2; ++ i) {
        edgeSlns.push_back( make_tuple(minElement+((i*epsilon)/2.0), BigInt(0)) );
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
        vector<size_t, StackAllocator<size_t>> inds;
        std::tie(cost, inds) = pq.front();
        std::pop_heap( pq.begin(), pq.end(), ord );
        pq.pop_back();

        // Compute the number of ways we can obtain this solution
        BigInt numSlns(1);
        for ( size_t i = 0; i < inds.size(); ++i ) {
            if (printMe) {
                auto vert = H->vertex(tailNodes[i]);
                cerr << vstr( vert ) << "\n";
                cerr << "i = " << i << ", score = " << countDict[ tailNodes[i] ][ inds[i] ].score << ", count = " << countDict[ tailNodes[i] ][ inds[i] ].count << "\n";
            }

            numSlns *= countDict[ tailNodes[i] ][ inds[i] ].count;
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
    BigInt zero(0);
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
        [&]( const tuple<double, BigInt>& e1, const tuple<double, BigInt>& e2 ) -> bool {
            return get<0>(e1) < get<0>(e2);
        }
    );
    if(edgeSlns.size() == k+4) { for (size_t i = 0; i < 4; ++i) { edgeSlns.pop_back(); } }
    */
}

vector<double> computeAlphasDouble( const double &bscale, const vector<MultiOpt::ScoreCount> &slnVec, size_t k, const BigInt &total ) {
    if (slnVec.size() == 0) {
        return vector<double>();
    }
    vector<double> scores;
    scores.reserve(slnVec.size());
    double bestScore = slnVec.front().score;
    double worstScore = slnVec.back().score;
    if ( bestScore == worstScore && slnVec.size() > 2 ) {
        cerr << "bestScore (" << bestScore << ") == worstScore (" << worstScore << ")" << "\n";
        cerr << "=== slnVec ===\n";
        for (const auto & e : slnVec) {
            cerr << "score = " << e.score << ", count = " << e.count << "\n";
        }
        std::abort();
    }
    double diff = (worstScore == bestScore) ? 1.0 : worstScore - bestScore; // std::max(0.01, worstScore - bestScore);
    size_t N = slnVec.size();

    //double scale = estimateThermodynamicBeta( slnVec, bestScore ); // (6.9*k) / diff; // (1.25 * k) / diff;
    double scale = bscale / diff;
    //double scale = bscale;

    double sum(0.0);
    size_t i = 0;
    for (const auto & e : slnVec) {
        double a = std::exp( -std::abs( (bestScore - e.score) * scale) );
        scores.emplace_back( a ); 
        sum += scores.back();
        i += 1;
    }

    auto invSum = 1.0 / sum;
    vector<double> alphas;
    alphas.reserve(slnVec.size());
    for (const auto & s : scores) {
        alphas.emplace_back( s * invSum );
    }

    return alphas;
}

template <typename CostClassT>
vector<CostClassT> computeKBest(const size_t &vid,
                                const size_t &k,
                                vector< vector<CostClassT> > &tkd,
                                unique_ptr<ForwardHypergraph> &H) {

    typedef vector<size_t, StackAllocator<size_t>> IndexVector;
    typedef vector<size_t, StackAllocator<size_t>> SizeVector;
    using MultiOpt::EdgeDerivation;

    // Dictionary that holds, for each incoming edge, the
    // number of score classes for each tail node
    unordered_map<size_t, IndexVector> sizeDict;

    // Priority queue of derivations for the given vertex
    using boost::heap::skew_heap;
    //QueueCmp<edvsT> ord;
    //skew_heap<edvsT, boost::heap::compare<QueueCmp<edvsT>>> vpq(ord);
    CountedDerivCmp<EdgeDerivation> ord;
    skew_heap<EdgeDerivation, boost::heap::compare<CountedDerivCmp<EdgeDerivation>>> vpq(ord);

    // Will hold the top-k cost classes for this solution
    vector<CostClassT> cc;
    cc.reserve(k);

    
    for ( auto & eid : H->incident(vid) ) {
        // The edge, backpointer array and cost of the derivation
        auto edge = H->edge(eid);
        IndexVector bp(edge.tail().size(), 0);
        auto cost = round3(getCost(tkd, bp, eid, H));

        // Again, check Boost emplace bug [Ticket #8195]
        vpq.emplace( cost, eid, bp );

        // Fill in the size array for this edge
        SizeVector sizes(edge.tail().size(), 0);
        size_t i = 0;
        for ( auto & tn : edge.tail() ) {
            sizes[i] = tkd[tn].size(); ++i;
        }
        sizeDict[eid] = std::move(sizes);
    }

    // Compute the cost of a derivation
    std::function< double( const size_t &, const IndexVector& ) > computeScore = [&] (const size_t & eid,
    const IndexVector &inds ) -> double {
        return round3( getCost(tkd, inds, eid, H) );
    };

    // Exact score classes
    double epsilon = 0.25;//0.1;

    // While there are still derivations left in the heap,
    // and we don't yet have the required number of solutions
    while ( not vpq.empty() and cc.size() <= k ) {

        // Get the top edge derivation
        const EdgeDerivation& ederiv = vpq.top();
        
        // Determine the number of solutions it yields
        auto count = getCount( tkd, ederiv.backPointers, ederiv.edgeID, H );
        // Construct the counted derivation
        CountedDerivation cderiv( ederiv.cost, ederiv.edgeID, ederiv.backPointers, count );
        // Now we no longer need ederiv, so pop it from the queue
        vpq.pop();

        // If the list of cost classes is empty, or if this
        // derivation belongs in a new cost class, then
        if ( cc.size() == 0 or (fabs(cderiv.cost - cc.back().cost())) > epsilon ) {
            // Create the new cost class and append the counted derivation
            CostClassT cclass(cderiv.cost);
            cclass.appendToCount(cderiv);
            cc.emplace_back( std::move(cclass) );
        } else { // we already have a solution of this score
            // Append the counted derivation to the last cost class
            cc.back().appendToCount(cderiv);
        }

        // Append the successors of this derivation to the heap
        Utils::appendNextWithEdge( cderiv.edge, cderiv.bp, sizeDict[cderiv.edge], vpq, computeScore );
    }

    // The cost class, the whole cost class and nothing
    // but the cost class
    if ( cc.size() == k + 1 ) {
        cc.resize(k);
    }
    return cc;
}

template <typename CostClassT>
bool viterbiCountNew( 
    unique_ptr<ForwardHypergraph> &H, //<! The problem hypergraph
    TreePtrT &t,                      //<! Phylogeny of the proteins
    TreeInfo &ti,                     //<! Extra information about the tree
    double penalty,                   //<! Penalty for creating edges that disagree with branch lengths
    const vector<size_t> &order,      //<! Topological order of hypergraph nodes
    slnDictT &slnDict,                //<! Holds solutions
    countDictT &countDict,            //<! Holds the number of ways of getting each vertex & cost 
    const size_t &k,
    const string &outputName, 
    const vector<FlipKey> &outputKeys, 
    const double &beta 
    ) {

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
            CountedDerivation cderiv( slnDict[vit][0].cost, std::numeric_limits<size_t>::max(), vector<size_t, CustomAllocator<size_t>>(), BigInt(1) );
            tkd[vit].back().appendToCount( cderiv );
            countDict[vit].emplace_back( slnDict[vit][0].cost, 1 );
            costsum += slnDict[vit][0].cost;
        }
    }

    cerr << "COSTSUM = " << costsum << "\n";
    cerr << "ORDER SIZE = " << order.size() << "\n";
    cerr << "SIZE OF COUNT DICT = " << countDict.size() << "\n";
    

    typedef size_t edgeIdT;
    auto N = order.size();

    ez::ezETAProgressBar eta(order.size()); eta.start();

    // For each vertex, in topological order (up)
    for ( auto vit = order.begin(); vit != order.end(); ++vit, ++eta ) {

        if ( H->incident(*vit).size() != 0 ) {
            auto vert = H->vertex(*vit);
            tkd[*vit] = std::move( computeKBest(*vit, k, tkd, H) );
        }

        // Find the minimum cost edge
        Google<Derivation::flipT>::Set fs;
        fs.set_empty_key( Derivation::flipT(-1, -1, "") );
        slnDict[*vit][0] = Derivation( tkd[*vit].front().cost(), 0, vector<size_t>(), fs);

    } // loop over verts

    typedef std::vector<double> ProbMap;
    typedef std::vector<size_t> NumTailMap;
    ProbMap probMap(order.size(), 0.0);
    ProbMap outProbMap(order.size(), 0.0);
    
    // NumTailMap ntmap(order.size(),0);
    // size_t invalidIdx = std::numeric_limits<size_t>::max();
    // typedef Google< size_t, double >::Map ProbMap;
    // ProbMap probMap;
    // ProbMap outProbMap;
    // probMap.set_empty_key( invalidIdx );
    // outProbMap.set_empty_key( invalidIdx );

    // Map from a vertex to the maximum cost class that is
    // considered in deriving any solution that is actually used.
    vector<size_t> maxCostClass(maxID + 1, 0);
    FlipKey rootKey( t->getRootId(), t->getRootId(), FlipState::none, MustConnect::none );
    auto rootInd = H->index(rootKey);

    // We always consider all k classes for the root. NOTE: There
    // may be less than k classes if we have enumerated all
    // of them.
    maxCostClass[rootInd] = std::min(k, tkd[rootInd].size());

    // The root gets a probability of 1
    probMap[rootInd] = 1.0;
    // Testing
    // ntmap[rootInd] = 1;

    size_t ctr = 0;
    size_t tot = order.size();
    auto rootFlip = flipBoth(rootKey);
    H->addVertex( rootFlip );
    auto rootIdNoFlip = H->index(rootKey);
    auto rootIdFlip = H->index(rootFlip);

    cerr << "Down phase\n";
    eta.reset(order.size()); eta.start();

    /******
    * Testing stuff
    */
    /*
    boost::dynamic_bitset<> normed(order.size()); // all 0's by default
    boost::dynamic_bitset<> normedOut(order.size()); // all 0's by default

    class ProbTransfer {
    private:
        size_t ind1;
        size_t ind2;
        std::vector<double> transferedMass;
    public:
        ProbTransfer() : ind1(0), ind2(0), transferedMass(std::vector<double>(4, 0.0)) {}
        ProbTransfer( size_t a, size_t b ) : ind1(a), ind2(b), transferedMass(std::vector<double>(4, 0.0)) {}
        ProbTransfer( size_t a, size_t b, std::vector<double>& v ) : ind1(a), ind2(b), transferedMass(v) {}

        void transferMass( size_t from, size_t to, double amt ) {
            if ( from == ind1 and to == ind1 ) {
                transferedMass[0] += amt;
            } else if ( from == ind1 and to == ind2 ) {
                transferedMass[1] += amt;
            } else if ( from == ind2 and to == ind1 ) {
                transferedMass[2] += amt;
            } else if ( from == ind2 and to == ind2 ) {
                transferedMass[3] += amt;
            }
        }

        void print( const std::string& k1, const std::string& k2 ) {
            std::cerr << "Prob transfer:\n" 
                      << k1 << " => " << k1 << " = " << transferedMass[0] << "\n"
                      << k1 << " => " << k2 << " = " << transferedMass[1] << "\n"
                      << k2 << " => " << k1 << " = " << transferedMass[2] << "\n"
                      << k2 << " => " << k2 << " = " << transferedMass[3] << "\n";
        }
    };

    std::unordered_map< std::tuple<size_t, size_t>, ProbTransfer > ptsMap;
    */


    /******
    * End testing stuff
    */

    // Compute the probabilities (down)
    // Loop over all vertices in *reverse* topological order
    for ( auto vit = order.rbegin(); vit != order.rend(); ++vit, ++ctr, ++eta ) {

        // The current vertex and it's probability
        auto key = H->vertex(*vit);
        auto parentProb = probMap[*vit];
        auto complementVert = flipBoth(key);
        auto complementInd = H->index(complementVert);

        /*
        // normalize the probabilities
        auto totalProb = probMap[*vit] + probMap[complementInd];
        
        if ( normed[*vit] == 0 ) {
           
           // re-adjust the weights
           if (totalProb > 0.0) {     
               probMap[*vit] /= totalProb;
               probMap[complementInd] /= totalProb;
           } else {
               probMap[*vit] = 0.5; outProbMap[*vit] = 0.5;
               probMap[complementInd] = 0.5; outProbMap[complementInd] = 0.5;
           }

           normed[*vit] = 1;
           normed[complementInd] = 1;
        } else {
            if( std::abs( totalProb - 1.0 ) >= 1e-10 ){
                std::cerr << "EDGE NOT NORMED; totalProb = " << totalProb << "\n"; 
                std::exit(0);
            }
        }
        
        parentProb = probMap[*vit];
        */
       
        // The total # of derivations of this vertex (over all considered cost classes)
        BigInt total(0);
        vector< MultiOpt::ScoreCount > cd;

        auto maxDeriv = maxCostClass[*vit];

        for ( size_t i = 0; i < maxDeriv; ++i ) {
            total += tkd[*vit][i].total();
            for ( auto & e : tkd[*vit][i].usedEdges() ) {
                auto& frontier = tkd[*vit][i].getEdgeFrontier(e);
                auto& tail = H->edge(e).tail();

                for ( size_t j : boost::irange(size_t{0}, tail.size()) ) {
                    auto tn = tail[j];
                    // eager
                    maxCostClass[tn] = std::max( frontier[j] + 1, maxCostClass[tn] );

                    // lazy
                    // maxCostClass[tn] = std::max( frontier[j][0]+1, maxCostClass[tn] );
                }
            }
            cd.emplace_back( tkd[*vit][i].cost(), tkd[*vit][i].total() );
        }

        // Compute the weights for each of the score classes
        // considered at this node
        auto alphas = computeAlphasDouble( beta, cd, maxCostClass[*vit], total );
        /** testing **/
        /*
        double asum = 0.0; for( auto a : alphas ) { asum += a; }
        if ( maxDeriv > 0 and std::abs(asum - 1.0) >= 1e-3 ) { 
            std::cerr << "ALPHA sum  = " << asum << " != 1!\n"; std::exit(1); 
        }
        auto totalOutProb = 0.0;
        */
       
        // for each top-k score class
        for ( size_t i = 0; i < maxDeriv; ++i ) {
            // The i-th cost class for this node
            auto& cc = tkd[*vit][i];
            auto scoreClassWeight = std::log(alphas[i]);
            // testing
            //double tprob = 0.0;

            // for all incoming edges contributing to this cost class
            for ( const auto & e : cc.usedEdges() ) {

                // The conditional probability of taking edge 'e'
                // given than we're deriving vertex *vit at score
                // the given cost
                auto condProb = cc.edgeProb(e);
                //tprob += std::exp(condProb);
                auto tail = H->getTail(e);

                // What is the action along edge 'e' (i.e. how is
                // the state of the interaction function different
                // between *vit and the tail nodes of e)
                auto ft = flipType( key, H->vertex(tail.front()) );
                auto outInd = (ft == "n") ? *vit : complementInd;

                // Accumulate the probability due to the current
                // cost class at outInd
                auto etprob = ( parentProb * (scoreClassWeight * condProb));

                // BEGIN Testing
                //double logPProb = std::log(parentProb);
                //double etprob = std::exp(logPProb + scoreClassWeight + condProb);
                //totalOutProb += etprob;
                // END Testing
                
                outProbMap[outInd] += etprob;

                /** Testing **/
                /*
                auto ptsKey = (*vit < complementInd) ? 
                              make_tuple(*vit, complementInd) : make_tuple(complementInd, *vit);

                auto pt = ptsMap.find(ptsKey);
                if ( pt == ptsMap.end() ) {
                    //std::cerr << "HERE\n";
                    //std::cerr << "keys " << get<0>(ptsKey) << ", " << get<1>(ptsKey) << "\n";
                    ptsMap.emplace(ptsKey, ProbTransfer(get<0>(ptsKey), get<1>(ptsKey)));
                    pt = ptsMap.find(ptsKey);
                }
                pt->second.transferMass( get<0>(ptsKey), get<1>(ptsKey), etprob );
                */
                /** End Testing **/

                // For each tail vertex of this hyperarc,
                // accumulate probability mass at this vertex
                for ( const auto & tind : tail ) {
                    probMap[tind] += etprob;
                    //probMap[tind] += tnprob;
                }
            } // end of loop over cost classes
            /** Testing **/
            /*
            if( not H->isLeaf(*vit) and std::abs(tprob - 1.0) >= 1e-5 ) { 
                std::cerr << "Conditional prob = " << tprob << "; NOT 1"; 
                
            }
            */
            /** End Testing**/
        } // end of loop over score classes
        /** Testing **/
        /*
        if( maxDeriv > 0 and (not H->isLeaf(*vit)) and std::abs(totalOutProb - parentProb) > 1e-5 ) { 
            std::cerr << "Total out prob (= " << totalOutProb << ") != "
                      << "Parent prob (= " << parentProb << ")\n"; 
        }
       
       
        normedOut[*vit] = 1;
        if ( (not H->isLeaf(*vit)) and normedOut[complementInd] == 1 ) {
            auto inSum = probMap[*vit] + probMap[complementInd];
            auto outSum = outProbMap[*vit] + outProbMap[complementInd];

            // The amount of mass coming *out* of the pair (x, complement(x))
            // is the same as the total mass going into it.
            if ( std::abs(outSum - inSum) > 1e-3 ) {
                std::cerr << "OUT SUM ERROR; Sum is " << outSum << "\n";
                std::cerr << "V = " << vstr(key) << " : Prob In " 
                          << probMap[*vit] << ", Prob Out " << outProbMap[*vit] << "\n";
                std::cerr << "compV = " << vstr(complementVert) << " : Prob In " 
                          << probMap[complementInd] << ", Prob Out " << outProbMap[complementInd] << "\n";

                std::cerr << "INFO\n";
                std::cerr << "=================\n";
                auto ptsKey = (*vit < complementInd) ? 
                              make_tuple(*vit, complementInd) : make_tuple(complementInd, *vit);
                ptsMap[ptsKey].print( vstr(H->vertex(get<0>(ptsKey))), vstr(H->vertex(get<1>(ptsKey))) );
                std::cerr << "=================\n";
            }
        }
        /** End Testing**/
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

    std::cerr << "Order size was " << order.size() << " probMap size was " << probMap.size();

    // for ( auto kv : probMap ) {
    //     size_t vid = kv.first;
    for ( size_t vid : boost::irange(size_t{0}, order.size()) ) {
        
        auto key = H->vertex(vid);
        if ( H->incident(vid).size() == 0 && vid != rootIdFlip && vid != rootIdNoFlip ) {
            if ( outProbMap[vid] != probMap[vid] ) { //&& outProbMap[vid] > probMap[vid] ) {
                LOG_WARN(log) << "node " << vstr(key) << ", inProbMap has " << probMap[vid] << ", outProbMap has " << outProbMap[vid] << "\n";
            }
            outProbMap[vid] = probMap[vid];
        }
    }


    for ( size_t vid : boost::irange(size_t{0}, order.size()) ) {
        auto& key = H->vertex(vid);
        auto fkey = flipBoth(key);
        auto oid = H->index(fkey);
        
        bool normalize = true;
        if ( normalize ) {
            auto psum = outProbMap[vid] + outProbMap[oid];        
            if ( psum > 0.0 )  {
                auto norm = 1.0 / psum;
                outProbMap[vid] *= norm;
                outProbMap[oid] *= norm;
            } else {
                outProbMap[vid] = outProbMap[oid] = 0.5;
            }

        }
        
        auto approxInProb = probMap[vid];
        auto approxProb = outProbMap[vid];

        bool writeOut = true;
        if ( restrictOutput ) {
            writeOut = ( std::find( outputKeys.begin(), outputKeys.end(), key ) != outputKeys.end() );
        }

        if ( approxProb > 0.0 && writeOut ) {
            auto fs = flipStrMap[key.state()];//.find(key.getDirTuple())->second;
            if ( key.state() != FlipState::none ) {
                output << t->getNodeName(key.u()) << "\t" << t->getNodeName(key.v())
                       << "\t" << fs << "\t" << approxProb << "\n";
            }
        }
    }
    output.close();
    LOG_INFO(log) << "leaving viterbiCountNew\n";
    std::cerr << "\n";

    //boost::singleton_pool<boost::pool_allocator_tag, sizeof(size_t)>::release_memory();
    return true;
}

vector<cl_RA> computeAlphas( const vector<tuple<double, BigInt>> &slnVec ) {
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

vector<double> getAlphas( const vector< tuple<double, BigInt> > &slnVec, const BigInt &numPaths ) {
    typedef tuple<double, BigInt> elemT;

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
    mpq_rational a;
    
    for ( const auto & e : slnVec) {
        //mpq_set_num(a.backend().data(), get<1>(e).backend().data() );//mpq_rational(get<1>(e), numPaths).convert_to<double>()
        //mpq_set_den(a.backend().data(), numPaths.backend().data() );
        //invScores.push_back( (alpha * a.convert_to<double>()) + ((1.0 - alpha) * std::exp( bestScore - get<0>(e) * scale ))  );

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

FlipKey canonicalKey( const FlipKey&k ){
    return FlipKey(k.u(), k.v(), k.state(), MustConnect::none);
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



double computeVertexProbability(const size_t &vid,
                                const TreePtrT& t,
                                const std::vector<double> &probs,
                                const boost::dynamic_bitset<> &normed,
                                unique_ptr<ForwardHypergraph> &H,
                                Model &model) {

    auto printWithNames = [&] (size_t ind) -> void {
        auto k = H->vertex(ind);
        auto dirStr = (k.f() and k.r()) ? " <--> " : " X ";
        std::cerr << "[ " << t->getNodeName(k.u()) << ", " << t->getNodeName(k.v()) << 
            " : " << dirStr << "] ";
    };

    /**
    * Determines if the transition from the key parent to the key child is valid.
    * Any transition in which the state of the interaction is not altered (not flipped) between the 
    * parent and child is valid.  If the state of the interaction is flipped, the transition is only
    * valid if the interacting nodes ( (u,v) in the parent ) exist in the same species.
    */
    auto isValidTransition = [&]( const FlipKey& parent, const FlipKey& child ) -> bool {
        if ( (parent.state() != child.state()) ) {
            auto us = dynamic_cast<bpp::BppString*>( t->getNodeProperty(parent.u(),"S") )->toSTL();
            auto vs = dynamic_cast<bpp::BppString*>( t->getNodeProperty(parent.v(),"S") )->toSTL();
            return us == vs;
        } 
        return true;        
    };

    auto incomingEdges = H->incident(vid);
    // We already know the probability of all of the leaves
    if ( incomingEdges.size() == 0 ) { return probs[vid]; }

    auto parentVertex = H->vertex(vid);
    double prob = 0.0;
    size_t k = 0;

    // For each edge incident to this vertex
    for (auto & eid : incomingEdges) {
        auto edge = H->edge(eid);

        // probability of the tail vertices
        double tprob = 1.0;
        size_t i = 0;

        // For every tail vertex of this hyperedge
        for (auto tid : edge.tail()) {
            auto childVertex = H->vertex(tid);
            //assert( normed[tid] );        
            tprob *= probs[tid] * model.transitionProbability( childVertex, parentVertex );
        }

        prob += tprob;
    }

    return prob;

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
    boost::dynamic_bitset<> normed( maxID + 1 );

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
            normed[vit] = 1;
        }
    }

    typedef size_t edgeIdT;
    auto N = order.size();
    size_t ctr = 0;
    ProgressDisplay showProgress(order.size());

    auto allNodesInvolving = [&] ( int u, int v ) -> std::vector<size_t> {
        std::vector<FlipState> fr = { FlipState::none, FlipState::both };
        std::vector<MustConnect> cucv = { MustConnect::none, MustConnect::left, MustConnect::right, MustConnect::both };
        std::vector< size_t > keys;

        for ( auto& dir : fr ) { // for all flip states
            for ( auto& cstate : cucv ) { // for all connect states
                FlipKey k(u, v, dir, cstate);
                if ( H->contains( k ) ) {
                    keys.push_back( H->index(k) );
                }
            }
        }
        return keys;
    };

    // For each vertex, in topological order (up)
    for ( auto vit = order.begin(); vit != order.end(); ++vit, ++ctr, ++showProgress ) {

        /*
         for ( auto e : H->incident(*vit) ) {
            auto edge = H->edge(e);
            for ( auto v : edge.tail() ) {
               auto key = H->vertex(v);
               
               auto opKey = flipBoth(key);
               auto opInd = H->index(opKey);
               
               if ( !H->isLeaf(v) ) {

                auto states = {v, opInd};
                double probSum = 0.0;
                for ( auto s : states ) { probSum += probs[s]; }                
                if ( ! (probSum > 0) ) {
                    for ( auto s : states ) {
                        if ( H->incident( opInd ).size() > 0 ) {
                            std::cerr << "P("; printWithNames( s ); std::cerr << ") = " << probs[s] << "\n";
                        }
                    }
                    //std::abort();
                 }
                 if ( probSum > 0 ) {
                    for ( auto s : states ){ 
                        probs[s] /= probSum;
                    }
                    //probs[v] /= probSum;
                    //probs[opInd] /= probSum;
                    //assert( probs[v] + probs[opInd] > 0.99 && probs[v] + probs[opInd] < 1.01 );
                    auto tmpSum = 0.0;
                    for ( auto s : states ) { tmpSum += probs[s]; }
                    assert( tmpSum > 0.99 && tmpSum < 1.01 );
                    
                 }
                 //normed[v] = true; normed[opInd] = true;
                 for ( auto s : states ) { normed[s] = 1; }
               }
           }
        } 
        */
       
        auto vert = H->vertex(*vit);
        auto probVert = computeVertexProbability(*vit, t, probs, normed, H, model);
        probs[ *vit ] = probVert;

    } // loop over verts

    typedef Google< size_t, double >::Map probMapT;
    size_t invalidIdx = std::numeric_limits<size_t>::max();
    probMapT probMap;
    probMapT outProbMap;
    probMap.set_empty_key( invalidIdx );
    outProbMap.set_empty_key( invalidIdx );

    // Map from a vertex to the maximum cost class that is
    // considered in deriving any solution that is actually used.
    FlipKey rootKey( t->getRootId(), t->getRootId(), FlipState::none );
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
            /*
            std::cerr << "key ";
            printWithNames(*vit);
            std::cerr << " is not normed!!\n";
            std::abort();
            */
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
            auto fs = flipStrMap[key.state()];//.find(key.getDirTuple())->second;
            if ( fs != "n" ) {
                output << t->getNodeName(key.u()) << "\t" << t->getNodeName(key.v())
                       << "\t" << fs << "\t" << approxProb << "\n";
            }
        }
    }
    output.close();
    std::cerr << "\n";

}


/* =============== BEGIN DEPRECATED ================== */

double estimateThermodynamicBeta( const vector<ScoreCount> &slnVecIn, const double &emin ) {
    if (slnVecIn.size() == 1) {
        return 0.0;
    }
    vector< tuple<double, mpreal> > slnVec;
    slnVec.reserve(slnVecIn.size());
    std::stringstream ss;
    for ( const auto & e : slnVecIn ) {
        ss << e.count;
        mpreal n = ss.str();
        ss.clear();
        slnVec.emplace_back( std::forward_as_tuple(e.score, n) );
    }
    double KB = 1.0;
    double beta = 0.0;
    size_t n = slnVec.size();
    double scale = 1.0 / ( n - 1 );
    
    for ( auto i = slnVec.begin(); (i + 1) != slnVec.end(); ++i ) {
        auto j = i + 1;
        auto lnI = KB * mpfr::log( get<1>(*i) );
        auto eI = get<0>(*i) - emin;
        auto lnJ = KB * mpfr::log( get<1>(*j) );
        auto eJ = get<0>(*j) - emin;
        double denom = eJ - eI;

        if ( denom > 0.0 ) {
            beta += scale * ((lnJ - lnI).toDouble() / denom);
        }
    }
    
    //beta /= ctr;
    std::cerr << " ===== beta = " << beta << "; num score classes = " << slnVecIn.size() << " ===== \n";
    return beta;
}


void insideOutside( unique_ptr<ForwardHypergraph> &H, TreePtrT &t, TreeInfo &ti, double penalty, const vector<size_t> &order, slnDictT &slnDict, countDictT &countDict, const string &outputName ) {

    // We'll use vectors for these since the nodes and vectors have a well
    // defined ordering.
    vector< tuple<double, BigInt> > insideScores( H->order() );
    vector< tuple<double, BigInt> > edgeWeightMap( H->size() );

    // Each leaf has a single solution which is, by definition, of optimal cost
    for ( const auto & vit : order ) {
        if ( H->incident(vit).size() == 0 ) {
            insideScores[vit] = make_tuple(slnDict[vit][0].cost, BigInt(1));
        } else {
            insideScores[vit] = make_tuple(0.0, BigInt(0));
        }
    }

    for ( size_t i = 0; i < H->size(); ++i ) {
        edgeWeightMap[i] = make_tuple(0.0, BigInt(0));
    }

    size_t ctr = 0;
    size_t tot = H->order();
    // For each vertex, in topological order (inside)
    for ( auto vit = order.begin(); vit != order.end(); ++vit, ++ctr ) {
        cerr << "\r\rprocessing node " << ctr << "/" << tot;
        vector< tuple<double, BigInt> > scores;
        scores.reserve( H->incident(*vit).size() );
        if ( H->incident(*vit).size() > 0 ) { // For every non-leaf
            BigInt totalPaths(1);
            // loop over all incoming edges
            for ( const auto & e : H->incident(*vit) ) {

                auto tail = H->getTail(e);
                vector< tuple<double, BigInt> > tailScores;
                tailScores.reserve( tail.size() );

                BigInt numPaths(1);
                BigInt summedNumPaths(0);
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

    FlipKey rootKey( t->getRootId(), t->getRootId(), FlipState::none );
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
        auto fs = flipStrMap[key.state()];//.find(key.getDirTuple())->second;
        if ( approxProb > 0.0 && fs != "n" ) {
            output << t->getNodeName(key.u()) << "\t" << t->getNodeName(key.v())
                   << "\t" << fs << "\t" << approxProb << "\n";

        }
    }
    output.close();
    std::cerr << "\n";



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
            countDict[vit].emplace_back( slnDict[vit][0].cost, 1 );
            costsum += slnDict[vit][0].cost;
        }
    }


    auto cdname = "CDICT.txt";
    std::fstream cdout( cdname, std::fstream::out | std::fstream::trunc );
    cdout << countDict.size() << "\n";
    /*
    for ( auto ele : countDict ) {
        size_t vit = ele.first;
        auto vert = H->vertex(vit);
        cdout << vstr(vert) << "\t" << ele.second.size() << "\n";
        for ( auto tup : ele.second ) {
            cdout << std::fixed << std::setprecision(18) << get<0>(tup) << "\t" << get<1>(tup) << "\n";
        }
    }
    */

    cerr << "COSTSUM = " << costsum << "\n";
    cerr << "ORDER SIZE = " << order.size() << "\n";
    cerr << "SIZE OF COUNT DICT = " << countDict.size() << "\n";
    typedef size_t edgeIdT;
    // For each edge, count the number of solutions having each score
    unordered_map< edgeIdT, unordered_map< double, BigInt > > edgeCountMap;
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
        map< double, vector<tuple<size_t, BigInt> > > edgeCostMap;

        //cerr << "SLN FOR NODE " << vstr(vert) << "\n";
        out << vstr(vert);// << H->incident(*vit).size() << "\n";
        std::vector<BigInt> nedgesln;
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

            BigInt edgeSum(0);
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
                BigInt count;
                tie(score, count) = ent;
                auto edgeContrib = make_tuple(e, count);
                edgeCostMap[score].push_back( edgeContrib );
                edgeCountMap[ e ][ score ] = count;
            }
        }

        std::sort(nedgesln.begin(), nedgesln.end(),
        []( const BigInt & x, const BigInt & y) {
            return x < y;
        }
                 );

        for ( auto e : nedgesln ) {
            out << "\t" << e;
        }
        out << "\n";

        // If we traversed any edges
        if ( edgeCostMap.size() > 0 ) {
            typedef tuple<double, BigInt, size_t> edgeSlnT;
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
                BigInt numSln(0);

                // Update the information at the derived vertices
                for ( const auto & edgeCount : providingEdges ) {
                    size_t edgeInd;
                    BigInt count;
                    tie(edgeInd, count) = edgeCount;
                    usedEdges[*vit][score].insert( edgeInd );
                    // update the total number of solutions of the
                    // derived vertex
                    numSln += count;
                }
                // There are 'numSln' derivations yielding *vit at
                // a score of 'score'
                countDict[*vit].emplace_back( score, numSln );

                // Now that we have a total count for the derived
                // vertex, compute each edge's contribution
                for ( const auto & edgeCount : providingEdges ) {
                    size_t edgeInd;
                    BigInt count;
                    tie(edgeInd, count) = edgeCount;
                    //double edgeProb = mpq_rational( count / numSln ).convert_to<double>();
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

    FlipKey rootKey( t->getRootId(), t->getRootId(), FlipState::none );
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
        BigInt total(0);
        for ( size_t i = 0; i < countDict[*vit].size(); ++i ) {
            total += countDict[*vit][i].count;
        }

        // Compute the weights for each of the score classes
        // considered at this node
        auto alphas = computeAlphasDouble( beta, countDict[*vit], k, total );

        // for each top-k score class
        for ( size_t i = 0; i < countDict[*vit].size(); ++i ) {
            
            // The score and it's count
            auto pScoreCount = countDict[*vit][i];
            double pScore = pScoreCount.score;
            BigInt pCount = pScoreCount.count;
            

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
        cout << "score class " << i << " has " << sc.count << " solutions of score  " << std::fixed << std::setprecision(16) << sc.score << "\n";
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
        BigInt total(0);
        for ( size_t i = 0; i < countDict[*vit].size(); ++i ) { total += get<1>(countDict[*vit][i]); }

        // Compute the weights for each of the score classes
        // considered at this node
        auto alphas = computeAlphasDouble( countDict[*vit], k, total );


        // for each top-k score class
        for ( size_t i = 0; i < countDict[*vit].size(); ++i ) {
            double pScore; BigInt pCount;
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
            auto fs = flipStrMap[key.state()];//.find(key.getDirTuple())->second;
            if ( fs != "n" ) {
                output << t->getNodeName(key.u()) << "\t" << t->getNodeName(key.v())
                       << "\t" << fs << "\t" << approxProb << "\n";
            }
        }
    }
    output.close();
    std::cerr << "\n";

}



void viterbi( unique_ptr<ForwardHypergraph> &H, TreePtrT &t, TreeInfo &ti, double penalty, const vector<size_t> &order, slnDictT &slnDict ) {
    auto N = H->order();
    size_t ctr = 0;
    for ( auto vit = order.begin(); vit != order.end(); ++vit, ++ctr ) {
        cerr << "\r\rProcessed " << 100.0 * (static_cast<float>(ctr) / N) << "% of the vertices";
        auto vert = H->vertex(*vit);
        //cout << "Processing : [" << t->getNodeName(vert.u()) << ", " << t->getNodeName(vert.v()) << " : (" << vert.f() << vert.r() << ")]\t";
        double cval = std::numeric_limits<double>::max();
        //if ( slnDict.find(*vit) != slnDict.end() ) {
        if ( slnDict[*vit].size() == 0 ) {
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
            vector<size_t, CustomAllocator<size_t>> bp(tailNodes.size(), 0); // backpointer
            // array (all 0's)
            double cost; BigInt count;
            // ALLOCATOR uncommnet -- std::tie (cost, count) = getCostCount(derivs, bp, eit, H);
            vector<size_t, CustomAllocator<size_t>> mybp; mybp.reserve(tailNodes.size());
            std::copy(bp.begin(), bp.end(), mybp.begin());

            auto deriv = CountedDerivation(cost, eit, mybp, count);

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
            CountedDerivation cderiv( slnDict[vit][0].cost, std::numeric_limits<size_t>::max(), vector<size_t, CustomAllocator<size_t>>(), BigInt(1) );
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

    using boost::heap::pairing_heap;
    CountedDerivCmp<CountedDerivation> ord;
    cand[vid] = pairing_heap<CountedDerivation, boost::heap::compare<CountedDerivCmp<CountedDerivation>>>(ord);

    for ( auto e : H->incident(vid) ) {
        if (! derivs[vid].back().hasEdge(e) ) {
            auto edge = H->edge(e);
            auto tailNodes = edge.tail();
            auto tailSize = tailNodes.size();
            std::vector<size_t, StackAllocator<size_t>> bp(tailSize, 0);
            double cost;
            BigInt count;
            std::tie( cost, count ) = getCostCount(derivs, bp, e, H);
            vector<size_t, CustomAllocator<size_t>> mybp; mybp.reserve(tailNodes.size());
            std::copy(bp.begin(), bp.end(), mybp.begin());

            auto deriv = CountedDerivation(cost, e, mybp, count);
            cand[vid].push(deriv);
        }
    }
}

bool lazyNext(
    unique_ptr<ForwardHypergraph> &H,
    boost::heap::pairing_heap<CountedDerivation, boost::heap::compare<CountedDerivCmp<CountedDerivation>>> &localCandidates,
    size_t eind,
    const vector<size_t, CustomAllocator<size_t>> &j,
    size_t kp,
    DerivStoreT &derivs) {
    /*
    auto edge = H->edge(eind);
    auto tailNodes = edge.tail();
    auto headNode = H->vertex(edge.head());
    auto edgeCost = edge.weight();

    size_t i = 0;
    while ( i < tailNodes.size() ) {
        vector<size_t, CustomAllocator<size_t>> jp(j);
        jp[i] += 1;
        lazyKthBest(H, tailNodes[i], jp[i] + 1, kp, derivs);

        if ( jp[i] < derivs[tailNodes[i]].size() ) {
            double cost; BigInt count;
            std::tie(cost, count) = getCostCount(derivs, jp, eind, H);
            auto deriv = CountedDerivation(cost, eind, jp, count);
            localCandidates.push(deriv);
        }
        if (j[i] != 0) {
            return true;
        }
        i += 1;
    }*/
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

        vector< MultiOpt::ScoreCount > cd;
        // The total # of derivations of this vertex (over all
        // considered score classes)
        BigInt total(0);
        for ( auto & cc : derivs[vit] ) {
            cd.emplace_back( cc.cost(), cc.total() );
            total += cc.total();
        }

        // Compute the weights for each of the score classes
        // considered at this node
        auto alphas = computeAlphasDouble( beta, cd, derivs[vit].size(), total );

        // for each top-k score class
        for ( size_t i = 0; i < derivs[vit].size(); ++i ) {

            // The score and it's count
            double pScore(derivs[vit][i].cost());
            BigInt pCount(derivs[vit][i].total());

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
            auto fs = flipStrMap[key.state()];//.find(key.getDirTuple())->second;
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

/* =============== END DEPRECATED ================== */

}

using GraphUtils::undirectedGraphT;
using GraphUtils::directedGraphT;
using std::unique_ptr;

template void MultiOpt::MLLeafCostDict< undirectedGraphT  >( unique_ptr<ForwardHypergraph> & , Utils::Trees::TreePtrT & , undirectedGraphT & , bool , double , double , MultiOpt::slnDictT & );

template void MultiOpt::MLLeafCostDict< directedGraphT >( unique_ptr<ForwardHypergraph> & , Utils::Trees::TreePtrT & , directedGraphT & , bool , double , double , MultiOpt::slnDictT & );

template void MultiOpt::leafCostDict< undirectedGraphT  >( unique_ptr<ForwardHypergraph> & , Utils::Trees::TreePtrT & , Utils::TreeInfo& ti, undirectedGraphT & , bool , double , double , MultiOpt::slnDictT & );

template void MultiOpt::leafCostDict< directedGraphT >( unique_ptr<ForwardHypergraph> & , Utils::Trees::TreePtrT & , Utils::TreeInfo& ti, directedGraphT & , bool , double , double , MultiOpt::slnDictT & );

using MultiOpt::LazyCostClass;
using MultiOpt::EagerCostClass;

template std::tuple<double, MultiOpt::BigInt> MultiOpt::getCostCount( vector<vector<LazyCostClass>> &, const vector<size_t, StackAllocator<size_t>> &, const size_t &, unique_ptr<ForwardHypergraph> &);
template std::tuple<double, MultiOpt::BigInt> MultiOpt::getCostCount( vector<vector<EagerCostClass>> &, const vector<size_t, StackAllocator<size_t>> &, const size_t &, unique_ptr<ForwardHypergraph> &);

template double MultiOpt::getCost( vector<vector<LazyCostClass>> &, const vector<size_t, StackAllocator<size_t>> &, const size_t &, unique_ptr<ForwardHypergraph> &);
template double MultiOpt::getCost( vector<vector<EagerCostClass>> &, const vector<size_t, StackAllocator<size_t>> &, const size_t &, unique_ptr<ForwardHypergraph> &);

template MultiOpt::BigInt MultiOpt::getCount( vector<vector<LazyCostClass>> &, const vector<size_t, StackAllocator<size_t>> &, const size_t &, unique_ptr<ForwardHypergraph> &);
template MultiOpt::BigInt MultiOpt::getCount( vector<vector<EagerCostClass>> &, const vector<size_t, StackAllocator<size_t>> &, const size_t &, unique_ptr<ForwardHypergraph> &);

template bool MultiOpt::viterbiCountNew<CostClass<EdgeDerivInfoEager>>( unique_ptr<ForwardHypergraph> &H,
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
