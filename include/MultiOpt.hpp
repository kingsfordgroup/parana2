#ifndef MULTI_OPT_HPP
#define MULTI_OPT_HPP

#include <fstream>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include <Bpp/Phyl/Tree.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>
#include "utils.hpp"
#include "FHGraph.hpp"
#include <boost/assign/list_of.hpp>
#include <boost/functional/hash.hpp>
#include "Derivation.hpp"
#include <cln/rational.h>
#include <cln/integer.h>
#include <cln/float.h>
#include <cln/real.h>
#include <cln/rational_io.h>
#include <cln/integer_io.h>
#include <cln/float_io.h>
#include <limits>
#include <sstream>

/** Google's dense hash set and hash map **/
#include <google/dense_hash_set>
#include <google/dense_hash_map>
#include <google/sparse_hash_map>


namespace std {
    template<>
    class hash<tuple<int,int>> {
    public:
        std::size_t operator()(const tuple<int,int>& t) const {
            size_t seed = 0;
            boost::hash_combine(seed, get<0>(t));
            boost::hash_combine(seed, get<1>(t));
            return seed;
        }
    };
}

namespace std {
    template<>
    class hash<tuple<size_t,size_t>> {
    public:
        std::size_t operator()(const tuple<int,int>& t) const {
            size_t seed = 0;
            boost::hash_combine(seed, get<0>(t));
            boost::hash_combine(seed, get<1>(t));
            return seed;
        }
    };
}


namespace std {
    template<>
    class hash<tuple<bool,bool>> {
    public:
        std::size_t operator()(const tuple<bool,bool>& t) const {
            size_t seed = 0;
            boost::hash_combine(seed, get<0>(t));
            boost::hash_combine(seed, get<1>(t));
            return seed;
        }
    };
}

namespace MultiOpt {

    using cln::cl_I;
    using cln::cl_F;
    using cln::cl_RA;
    using cln::cl_R;
    using cln::cl_float;
    using google::dense_hash_set;
    using google::dense_hash_map;
    using google::sparse_hash_map;
    using std::map;
    using std::unordered_map;
    using std::unordered_set;
    using std::tuple;
    using std::tie;
    using std::make_tuple;
    using std::string;
    using std::get;
    using bpp::Tree;
    using Utils::TreeInfo;
    using Utils::differentExtantNetworks;
    using boost::assign::map_list_of;
    using std::unique_ptr;


    template< typename T1, typename T2=std::nullptr_t >
    struct Std{
        typedef unordered_set<T1, std::hash<T1>> Set;
        typedef unordered_map<T1, T2, std::hash<T1>> Map;
    };

    template< typename T1, typename T2=std::nullptr_t >
    struct Google{
        typedef dense_hash_set<T1, std::hash<T1>> Set;
        typedef dense_hash_map<T1, T2, std::hash<T1>> Map;
    };

    typedef tuple<bool,bool> IITup;
    typedef unordered_map<IITup,string> IITupStringMapT;
    typedef unordered_map<IITup, IITupStringMapT> FlipMapT;

    auto none = make_tuple(false,false); auto fw = make_tuple(true,false); auto rev = make_tuple(false,true); auto both = make_tuple(true,true);

    static unordered_map<tuple<bool,bool>, string> flipStrMap = {
        {none, "n"},
        {fw, "f"},
        {rev, "r"},
        {both, "b"}
    };

    typedef unordered_map< size_t, unordered_map<size_t, Derivation>> slnDictT;
    typedef Google< size_t, vector< tuple<double, cl_I> > >::Map countDictT;
    //typedef unordered_map<size_t, vector<Derivation>> slnDictT;
    //typedef Google<size_t, vector<Derivation>>::Map slnDictT;
    //typedef Google<size_t, Google<size_t,Derivation>::Map>::Map slnDictT;

    typedef tuple<bool,bool> flipTupleT;
    typedef tuple<double, string> costRepT;
    typedef unordered_map< flipTupleT, unordered_map<flipTupleT, costRepT> > costMapT;
    typedef unordered_map< flipTupleT, unordered_map<flipTupleT, std::function< costRepT (const double&, const double&) > > > costMapFunT;
    typedef unordered_map< flipTupleT, unordered_map<bool, costRepT> > selfCostMapT;
    typedef unordered_map< flipTupleT, unordered_map<bool, std::function<costRepT (const double&) > > > selfCostMapFunT;


    static FlipMapT flipDict = {
        {none, { {both, "b+"}, {fw, "f+"}, {rev, "r+"}, {none, "n"} } },
        {fw,   { {both, "r+"}, {rev, "f-r+"}, {none, "f-"}, {fw, "n"} } },
        {rev,  { {both, "f+"}, {rev, "n"}, {none, "r-"}, {fw, "f+r-"} } },
        {both, { {both, "n"}, {fw, "r-"}, {rev, "f-"}, {none, "b-"} } }
    };

    costMapT getCostDict ( double cc, double dc, bool directed ) {
        auto undirected = ! directed;
        auto none = make_tuple(false,false); auto fw = make_tuple(true,false);
        auto rev = make_tuple(false,true); auto both = make_tuple(true,true);

        costMapT costDict;
        costDict[ none ][ none ] = make_tuple( 0.0, "n" );
        costDict[ none ][ fw ] = make_tuple( cc, "f+" );
        costDict[ none ][ rev ] = make_tuple( cc, "r+" );
        costDict[ none ][ both ] = make_tuple( undirected ? cc : 2*cc , "b+" );

        costDict[ rev ][ none ] = make_tuple( dc, "r-" );
        costDict[ rev ][ fw ] = make_tuple( cc+dc, "f+r-" );
        costDict[ rev ][ rev ] = make_tuple( 0.0, "n" );
        costDict[ rev ][ both ] = make_tuple( cc , "f+" );

        costDict[ fw ][ none ] = make_tuple( dc, "f-" );
        costDict[ fw ][ fw ] = make_tuple( 0.0, "n" );
        costDict[ fw ][ rev ] = make_tuple( cc+dc, "f-r+" );
        costDict[ fw ][ both ] = make_tuple( cc , "r+" );

        costDict[ both ][ none ] = make_tuple( undirected ? dc : 2*dc, "b-" );
        costDict[ both ][ fw ] = make_tuple( dc, "r-" );
        costDict[ both ][ rev ] = make_tuple( dc, "f-" );
        costDict[ both ][ both ] = make_tuple( 0.0, "n");

        return costDict;
    }

    costMapFunT getCostFunDict ( double cc, double dc, bool directed ) {
        auto undirected = ! directed;
        auto none = make_tuple(false,false); auto fw = make_tuple(true,false);
        auto rev = make_tuple(false,true); auto both = make_tuple(true,true);

        costMapFunT costDict;
        costDict[ none ][ none ] = [=](const double& pf, const double& pr) { return make_tuple( 0.0, "n" ); };
        costDict[ none ][ fw ] = [=](const double& pf, const double& pr) { return make_tuple( pf*cc, "f+" ); };
        costDict[ none ][ rev ] = [=](const double& pf, const double& pr) { return make_tuple( pr*cc, "r+" ); };
        costDict[ none ][ both ] = [=](const double& pf, const double& pr) { return make_tuple( undirected ? pf*cc : (pf+pr)*cc, "b+" ); };


        costDict[ rev ][ none ] = [=](const double& pf, const double& pr) { return make_tuple( pr*dc, "r-" ); };
        costDict[ rev ][ fw ] = [=](const double& pf, const double& pr) { return make_tuple( pf*cc + pr*dc, "f+r-" ); };
        costDict[ rev ][ rev ] = [=](const double& pf, const double& pr) { return make_tuple( (1.0-pr)*dc, "n"); };
        costDict[ rev ][ both ] = [=](const double& pf, const double& pr) { return make_tuple( pf*cc + (1.0-pr)*dc , "f+"); };

        costDict[ fw ][ none ] = [=](const double& pf, const double& pr) { return make_tuple( pf*dc, "f-"); };
        costDict[ fw ][ fw ] = [=](const double& pf, const double& pr) { return make_tuple( (1.0-pf)*dc, "n"); };
        costDict[ fw ][ rev ] = [=](const double& pf, const double& pr) { return make_tuple( pf*dc + pr*cc, "f-r+"); };
        costDict[ fw ][ both ] = [=](const double& pf, const double& pr) { return make_tuple( (1.0-pf)*dc + pr*cc, "r+"); };

        costDict[ both ][ none ] = [=](const double& pf, const double& pr) { return make_tuple( undirected ? pf*dc : (pf+pr)*dc, "b-" ); };
        costDict[ both ][ fw ] = [=](const double& pf, const double& pr) { return make_tuple( (1.0-pf)*dc + pr*dc, "r-" ); };
        costDict[ both ][ rev ] = [=](const double& pf, const double& pr) { return make_tuple( pf*dc + (1.0-pr)*dc, "f-" ); };
        costDict[ both ][ both ] = [=](const double& pf, const double& pr) { return make_tuple( undirected ? (1.0-pf)*dc : (1.0-pf)*dc + (1.0-pr)*dc, "n" ); };

        return costDict;
    }

    selfCostMapT getSelfLoopCostDict ( double cc, double dc, bool directed ) {
        auto undirected = ! directed;
        auto none = make_tuple(false,false); auto fw = make_tuple(true,false);
        auto rev = make_tuple(false,true); auto both = make_tuple(true,true);

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
        auto none = make_tuple(false,false); auto fw = make_tuple(true,false);
        auto rev = make_tuple(false,true); auto both = make_tuple(true,true);

        selfCostMapFunT selfLoopCostDict;
        selfLoopCostDict[ none ][ true ] = [=] (const double& p) { return make_tuple( p*cc, "b+" ); };
        selfLoopCostDict[ none ][ false ] = [=] (const double& p) { return make_tuple( 0, "n" ); };

        selfLoopCostDict[ rev ][ true ] = [=] (const double& p) { return make_tuple( (1.0-p)*dc, "n" ); };
        selfLoopCostDict[ rev ][ false ] = [=] (const double& p) { return make_tuple( p*dc, "b-" ); };

        selfLoopCostDict[ fw ][ true ] = [=] (const double& p) { return make_tuple( (1.0-p)*dc, "n" ); };
        selfLoopCostDict[ fw ][ false ] = [=] (const double& p) { return make_tuple( p*dc, "b-" ); };

        selfLoopCostDict[ both ][ true ] = [=] (const double& p) { return make_tuple( (1.0-p)*dc, "n" ); };
        selfLoopCostDict[ both ][ false ] = [=] (const double& p) { return make_tuple( p*dc, "b-" ); };
        return selfLoopCostDict;
    }




    string flipType( const FlipKey& hvert, const FlipKey& tvert ){
        return flipDict[ make_tuple(hvert.f(), hvert.r())][ make_tuple(tvert.f(), tvert.r()) ];
    }

    template<typename GT>
    void projectToReversedGraph( unique_ptr<ForwardHypergraph>& H, GT& G ) {
        auto M = H->size();

        for ( size_t eid = 0; eid < M; ++eid ) {
            auto e = H->edge(eid);
            auto h = e.head();
            for ( auto t : e.tail() ) {
                add_edge(h, t, G);
            }
        }

    }

    void topologicalOrder( unique_ptr<ForwardHypergraph>& H, vector<size_t>& order ) {
        using boost::adjacency_list;
        using boost::vecS;
        using boost::directedS;

        typedef adjacency_list<vecS, vecS, directedS> graphT;
        graphT G;
        projectToReversedGraph( H, G );
        boost::topological_sort(G, std::back_inserter(order));
    }

    unique_ptr<ForwardHypergraph>  buildSolutionSpaceGraph( const Tree* t,
                                                            const TreeInfo& ti,
                                                            double cc,
                                                            double dc,
                                                            double penalty,
                                                            bool directed ) {
        Google<int>::Set leafSet;
        leafSet.set_empty_key(-1);
        for ( auto l : t->getLeavesId() ){ leafSet.insert(l); }
        auto isLeaf = [&]( const int& nid ) { return leafSet.find(nid) != leafSet.end(); };//t->isLeaf(nid); };
        auto isInternal = [&]( const int& nid ) { return ! isLeaf(nid); };
        auto isLost = [&]( int nid ){ return (t->getNodeName(nid)).find("LOST") != std::string::npos; };

        unique_ptr<ForwardHypergraph> slnSpaceGraph( new ForwardHypergraph() );

        costMapT costMap( getCostDict(cc,dc,directed) );
        selfCostMapT selfLoopCostMap( getSelfLoopCostDict(cc,dc,directed) );

        auto addIncomingHyperedge = [=,&costMap,&selfLoopCostMap,&slnSpaceGraph] (const FlipKey& k, const int& rnode, const int& onode ) {
            auto dirKey = k.getDirTuple();

            // Self loop
            if (k.arity() == 1) {
                auto selfNode = rnode;
                if (isInternal(selfNode)) {
                    // Get the children nodes

                    int LRN = t->getSonsId(selfNode)[0]; int RRN = t->getSonsId(selfNode)[1];
                    // This vertex has 2 incoming hyper edges
                    // 1 -- we don't flip the self-loop
                    auto noFlipLL = FlipKey( LRN, LRN, k.f(), k.r() );
                    auto noFlipRR = FlipKey( RRN, RRN, k.f(), k.r() );
                    auto noFlipLR = FlipKey( LRN, RRN, k.f(), k.r() );
                    // 2 -- we flip the self loop
                    auto dualFlipLL = flipBoth( noFlipLL ); auto dualFlipRR = flipBoth( noFlipRR ); auto dualFlipLR = flipBoth( noFlipLR );
                    vector<FlipKey> noFlipEdge = vector<FlipKey>();
                    vector<FlipKey> dualFlipEdge = vector<FlipKey>();

                    if ( !isLost(LRN) ) {
                        noFlipEdge.push_back( noFlipLL ); dualFlipEdge.push_back( dualFlipLL );
                    }
                    if ( !isLost(RRN) ) {
                        noFlipEdge.push_back( noFlipRR ); dualFlipEdge.push_back( dualFlipRR );
                    }


                    //vector<FlipKey> noFlipEdge = { noFlipLL, noFlipRR };
                    //vector<FlipKey> dualFlipEdge = { dualFlipLL, dualFlipRR };

                    if (! differentExtantNetworks(ti, LRN, RRN) && ! (isLost(LRN) || isLost(RRN)) ) {
                        noFlipEdge.push_back( noFlipLR );
                        dualFlipEdge.push_back( dualFlipLR );
                    }
                    if ( noFlipEdge.size() > 0 ) { slnSpaceGraph->addEdge( noFlipEdge, k, 0.0 ); }
                    if ( dualFlipEdge.size() > 0 ) {
                        slnSpaceGraph->addEdge( dualFlipEdge, k,
                                                get<0>(selfLoopCostMap[ make_tuple(k.f(), k.r()) ][ dualFlipLL.f() || dualFlipLL.r() ])  );
                    }
                }
            } else {
                if ( isInternal(rnode) ) {
                    //cout << t->getNodeName(rnode) << " is INTERNAL\n";
                    // Get the children nodes
                    int LRN = t->getSonsId(rnode)[0]; int RRN = t->getSonsId(rnode)[1];

                    double canonicalDerivCost = 0.0;
                    // if ( rnode < onode ) { canonicalDerivCost = 1e-50; }

                    // This vertex has 2 incoming hyper edges
                    // 1 -- we don't flip the self-loop
                    auto noFlipL = FlipKey( LRN, onode, k.f(), k.r() );
                    auto noFlipR = FlipKey( RRN, onode, k.f(), k.r() );
                    // 2 -- we flip the self loop
                    auto dualFlipL = flipBoth( noFlipL ); auto dualFlipR = flipBoth( noFlipR );

                    vector<FlipKey> noFlip; vector<FlipKey> dualFlip;

                    if ( !differentExtantNetworks(ti, LRN, onode) && !(isLost(LRN) || isLost(onode)) ){
                        noFlip.push_back( noFlipL ); dualFlip.push_back( dualFlipL );
                    }
                    if ( !differentExtantNetworks(ti, RRN, onode) && !(isLost(RRN) || isLost(onode)) ){
                        noFlip.push_back( noFlipR ); dualFlip.push_back( dualFlipR );
                    }

                    if ( noFlip.size() > 0 ) {
                        slnSpaceGraph->addEdge( noFlip, k, canonicalDerivCost);//0.0 );
                    }
                    if ( dualFlip.size() > 0 ) {
                        auto w = get<0>(costMap[ make_tuple(k.f(),k.r()) ][ make_tuple(dualFlipL.f(), dualFlipL.r()) ]);//directed ? 2.0 : 1.0;//2.0 ? directed : 1.0;
                        slnSpaceGraph->addEdge( dualFlip, k, w );
                    }

                    if ( directed ) {
                        auto fwFlipL = flipForward(noFlipL); auto fwFlipR = flipForward(noFlipR);
                        auto revFlipL = flipReverse(noFlipL); auto revFlipR = flipReverse(noFlipR);

                        vector<FlipKey> fwFlip; vector<FlipKey> revFlip;

                        if ( !differentExtantNetworks(ti, LRN, onode) ) {
                            fwFlip.push_back( fwFlipL ); revFlip.push_back( revFlipL );
                        }

                        if ( !differentExtantNetworks(ti, RRN, onode) ) {
                            fwFlip.push_back( fwFlipL ); revFlip.push_back( revFlipR );
                        }

                        if ( fwFlip.size() > 0 ) {
                            auto w = get<0>(costMap[ make_tuple(k.f(), k.r() ) ][ make_tuple(fwFlipL.f(), fwFlipL.r()) ]);
                            slnSpaceGraph->addEdge( fwFlip, k , w );
                        }
                        if ( revFlip.size() > 0 ) {
                            auto w = get<0>(costMap[ make_tuple(k.f(), k.r() ) ][ make_tuple(revFlipL.f(), revFlipL.r()) ]);
                            slnSpaceGraph->addEdge( revFlip, k, w );
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

        auto tbegin = nodes.cbegin(); auto tend = nodes.cend();
        vector<int>::const_iterator uit = tbegin; vector<int>::const_iterator vit = tbegin;


        for ( uit = tbegin; uit != tend; uit++ ) {
            int u = *uit;
            for ( vit = uit; vit != tend; vit++ ) {
                int v = *vit;

                if ( (u == v) ||
                     (!differentExtantNetworks(ti,u,v) &&
                      !(ti.inSubnodesOf(u,v) ||  ti.inSubnodesOf(v,u)) &&
                      !(isLost(u) || isLost(v))) ) {
                    //cout << "HERE\n";
                    slnSpaceGraph->addVertex( FlipKey( u, v, false, false ) );
                    if ( ! t->isRoot(v) ) {
                        slnSpaceGraph->addVertex( FlipKey( u, v, true, true ) );
                    }
                    if ( directed && u!=v ) {
                        slnSpaceGraph->addVertex( FlipKey( u, v, true, false ) );
                        slnSpaceGraph->addVertex( FlipKey( u, v, false, true ) );
                    }

                }
            }
        }

        auto N = slnSpaceGraph->order();
        cerr << "Hypergraph order = " << slnSpaceGraph->order() << "\n";
        for ( size_t i = 0; i < N; ++i ) {
            cerr << "\r\rProcessed " << i << "/" << N << " nodes";
            auto k = slnSpaceGraph->vertex(i);
            auto u = k.u(); auto v = k.v(); auto f = k.f(); auto r = k.r();
            addIncomingHyperedge( k, u, v );
            if ( k.arity() > 1 ) {
                addIncomingHyperedge( k, v, u );
            }
        }

        cerr << "\n";
        cerr << "Hypergraph size = " << slnSpaceGraph->size() << "\n";
        return slnSpaceGraph;
    }



    template< typename GT >
    void leafCostDict( unique_ptr<ForwardHypergraph>& H, Tree* T, GT& G, bool directed, double cc, double dc, slnDictT& slnDict ) {
        /*
          Given the duplication tree T, the root vertex rv, the extant graph G and
          the constraints, fill in the values for the leaf nodes of the hypergraph
        */
        typedef typename boost::graph_traits< GT >::edge_descriptor EdgeT;

        auto undirected = !directed;
        // Cost of going from the state encoded by a node to the state of a true graph edge

        auto none = make_tuple(false,false); auto fw = make_tuple(true,false);
        auto rev = make_tuple(false,true); auto both = make_tuple(true,true);

        auto costDict = getCostDict(cc,dc,directed);
        auto selfLoopCostDict = getSelfLoopCostDict(cc,dc,directed);
        auto costFunDict = getCostFunDict(1.0, 1.0, directed);
        auto selfLoopCostFunDict = getSelfLoopCostFunDict(1.0, 1.0, directed);

        auto N = H->order();
        auto M = H->size();

        cerr << "total # of hypernodes = " << N << "\n";
        cerr << "total # of hyperedges = " << M << "\n";

        // The list of all hypernodes with no descendants
        // We'll have either N^2 or (N^2)/2 leaf hypernodes (depending
        // on directedness)
        auto numExtantNodes = std::distance( vertices(G).first, vertices(G).second );
        auto numConn = ( numExtantNodes * numExtantNodes );
        if (undirected) { numConn /= 2; }
        vector<size_t> leafHypernodes; leafHypernodes.reserve( numConn);

        // For every hypernode, it's directed iff it has no incoming edges
        for( size_t i = 0; i < N; ++i) {
            auto elist = H->incident(i);
            if ( elist.size() == 0 ) { leafHypernodes.push_back(i); }
        }

        cerr << "# of leaf hypernodes = " << leafHypernodes.size() << "\n";

        typedef typename GT::vertex_descriptor NodeT;
        typedef unordered_set<NodeT> NodeSetT;

        // Does the graph have an edge u,v
        auto hasEdge = [&] ( size_t u, size_t v ) {
            auto av = boost::adjacent_vertices(u,G);
            for ( auto it = av.first; it != av.second; ++it ) {
                if ( *it == v ) { return true; }
            }
            return false;
        };

        // Is the node e contained in the set s?
        auto contains = [] ( const NodeSetT& s, NodeT e ) { return s.find(e) != s.end(); };

        NodeSetT extantNodes;

        Google<int>::Set leafIds;
        leafIds.set_empty_key(-1);

        for ( auto l : T->getLeavesId() ) { leafIds.insert(l); }

        auto vp = boost::vertices(G);
        for ( auto it = vp.first; it != vp.second; ++it ) {
            auto v = *it;
            auto idx = G[v].idx;
            // found this node's id in the set of extant vertices'
            if ( leafIds.find(idx) != leafIds.end() ) { extantNodes.insert(v); }
        }

        // Map from tree node ID to graph vertex ID
        unordered_map<int, NodeT> idToVertMap;
        for ( auto v = boost::vertices(G).first; v != boost::vertices(G).second; ++ v ){
            idToVertMap[ G[*v].idx ] = *v;
        }

        // For every leaf hypernode
        for( auto n : leafHypernodes ) {
            auto nd = H->vertex(n);

            auto nameU = T->getNodeName(nd.u());
            auto nameV = T->getNodeName(nd.v());

            // Check to see if u, v, or both have been lost
            auto endOfMap = idToVertMap.end();
            bool lostU = ( idToVertMap.find( nd.u() ) == endOfMap ) || ( ! contains(extantNodes, idToVertMap[nd.u()]) );
            bool lostV = ( idToVertMap.find( nd.v() ) == endOfMap ) || ( ! contains(extantNodes, idToVertMap[nd.v()]) );

            // The cost to / between lost nodes is always 0
            if ( lostU || lostV ) {
                auto lostCost = 0.0;
                vector<size_t> ev;
                Google<Derivation::flipT>::Set es;
                es.set_empty_key( make_tuple(-1,-1,""));
                slnDict[n] = { {0,Derivation(lostCost, n, ev, es)} };
            } else {
                // Otherwise, u and v both exist in the extant
                // network, so get the appropriate info
                auto u = idToVertMap[nd.u()];
                auto v = idToVertMap[nd.v()];
                auto f = nd.f(); auto r = nd.r();

                if (u != v) {
                    EdgeT fedge, redge;
                    bool d_f, d_r;
                    tie(fedge, d_f) = edge(u,v,G);
                    double w_f = d_f ? G[edge(u,v,G).first].weight : 0.0;
                    tie(redge, d_r) = edge(v,u,G);
                    double w_r = d_r ? G[edge(v,u,G).first].weight : 0.0;
                    if ( undirected ) {
                        assert( w_f == w_r );
                    }

                    //auto costFlip = costDict[ make_tuple(f,r) ][ make_tuple(d_f,d_r) ];
                    auto costFlip = costFunDict[ make_tuple(f,r) ][ make_tuple(d_f,d_r) ](w_f, w_r);

                    auto cost = get<0>(costFlip); auto flip = get<1>(costFlip);

                    Google<Derivation::flipT>::Set effectiveEdges;
                    effectiveEdges.set_empty_key( make_tuple(-1,-1,""));

                    if ( flip != "n" ) { effectiveEdges.insert( make_tuple(nd.u(),nd.v(),flip) ); }
                    vector<size_t> ev;
                    slnDict[n] = { {0,Derivation(cost, n, ev, effectiveEdges)} };

                } else {
                    EdgeT e; bool hasSelfLoop;
                    tie(e, hasSelfLoop) = edge(u,v,G);
                    double w_l = hasSelfLoop ? G[edge(u,v,G).first].weight : 0.0;

                    //auto costFlip = selfLoopCostDict[ make_tuple(f,r) ][ hasSelfLoop ];
                    auto costFlip = selfLoopCostFunDict[ make_tuple(f,r) ][ hasSelfLoop ]( w_l );

                    auto cost = get<0>(costFlip); auto flip = get<1>(costFlip);
                    Google<Derivation::flipT>::Set effectiveEdges;
                    effectiveEdges.set_empty_key( make_tuple(-1,-1,""));

                    if ( flip != "n" ) { effectiveEdges.insert( make_tuple(nd.u(),nd.v(),flip) ); }
                    vector<size_t> ev;
                    slnDict[n] = { {0,Derivation(cost, n, ev, effectiveEdges)} };

                } // ( u != v )
            } // ( lostU || lostV )
        } // loop over leaf hypernodes
    }

    /**
     *  Compute the penalty for this edge to exist based on difference
     *  between the existence intervals of the endpoints and the
     *  penalty factor.
     */
    template <typename T>
    double existencePenalty( TreeInfo& ti, const T& vert, const double& penalty ) {
        if (vert.f() || vert.r()) {
            if (ti.intervalDistance( vert.u(), vert.v() ) > 0.0) {
                return penalty * ti.intervalDistance( vert.u(), vert.v() );
            }
        }
        return 0.0;
    }

    typedef tuple<double, vector<size_t> > dvsT;

    class QueueCmp {
    public:
        bool operator() ( const dvsT& lhs, const dvsT& rhs ) {
            return get<0>(lhs) > get<0>(rhs);
        }
    };


    // 0 : 20 23 25
    // 1 : 15 16 34
    // 2 : 10 12 13
    // rank, #
    // (0, 3000)
    //
    vector< tuple<double, cl_I> > countEdgeSolutions( const double& ecost,
                                                      const vector<size_t>& tailNodes,
                                                      countDictT& countDict,
                                                      const size_t& k ) {
        // product pointers
        std::vector< size_t > prodElems(0, tailNodes.size());
        std::vector< size_t > elemSizes; elemSizes.reserve(tailNodes.size());
        double cost = ecost;
        for ( const auto& t : tailNodes ) {
            elemSizes.push_back( countDict[t].size() );
            cost += get<0>(countDict[t].front());
        }

        vector<dvsT> pq(1, make_tuple(cost, vector<size_t>(tailNodes.size(), 0)));
        QueueCmp ord;

        std::function< double( const vector<size_t>& ) > computeScore = [&] ( const vector<size_t>& inds ) -> double {
            size_t numNodes = tailNodes.size();
            double cost = ecost;
            for ( size_t i = 0; i < numNodes; ++i ) {
                cost += get<0>(countDict[ tailNodes[i] ][ inds[i] ]);
            }
            return cost;
        };

        typedef tuple<double, cl_I> ccT;
        vector< ccT > edgeSlns;
        double epsilon = 1e-5;
        size_t numSln = 0;

        while ( !pq.empty() && numSln < k ) {

            double cost;
            vector<size_t> inds;
            // Get the next best solution score from the top of the queue
            std:tie(cost, inds) = pq.front();
            // Compute the number of ways we can obtain this solution
            cl_I numSlns(1);
            for ( size_t i = 0; i < inds.size(); ++i ) { numSlns *= get<1>(countDict[ tailNodes[i] ][ inds[i] ]); }

            // put this solution into our # sln dict
            auto fp = std::find_if( edgeSlns.begin(), edgeSlns.end(), [=]( const ccT& cc ) ->bool { return (abs(cost - get<0>(cc)) <= epsilon);  } );
            if (fp == edgeSlns.end()) { // we didn't find a solution  of this score
                edgeSlns.push_back( make_tuple( cost, numSlns ) );
            } else { // we found a solution of this score
                get<1>(*fp) = get<1>(*fp) + numSlns;
            }

            Utils::appendNext( elemSizes, pq, ord, computeScore );
            numSln = edgeSlns.size();
        }
        return edgeSlns;
    }


    vector<cl_RA> computeAlphas( const vector<tuple<double, cl_I>>& slnVec ) {
        vector<cl_RA> scores; scores.reserve(slnVec.size());
        cl_RA invSum(0);
        cl_RA one(1);
        for (const auto& e : slnVec) {
            scores.push_back( one / static_cast<size_t>(get<0>(e)+1.0) );
            invSum += scores.back();
        }
        //cl_RA invSum = 1 / sum;
        // cerr << "invSum = " << invSum << "\n";
        vector<cl_RA> alphas; alphas.reserve(slnVec.size());
        for (const auto& s : scores) {
            alphas.push_back( s/invSum );
        }
        return alphas;
    }

    vector<double> computeAlphasDouble( const vector<tuple<double, cl_I>>& slnVec ) {
        vector<double> scores; scores.reserve(slnVec.size());
        double bestScore = get<0>(slnVec.front());
        double worstScore = get<0>(slnVec.back());
        double diff = std::max(1.0, worstScore - bestScore);
        double scale = 40. / diff;
        double invSum(0);
        //cl_RA one(1);
        for (const auto& e : slnVec) {
            scores.push_back( std::exp( (bestScore-get<0>(e)) * (scale) ) );/// (diff*diff) ) );//1.0 /( get<0>(e)+1.0 ) );
            invSum += scores.back();
        }

        vector<double> alphas; alphas.reserve(slnVec.size());
        for (const auto& s : scores) {
            alphas.push_back( s/invSum );
        }
        return alphas;
    }



    void viterbiCount( unique_ptr<ForwardHypergraph>& H, Tree* t, TreeInfo& ti, double penalty, const vector<size_t>& order, slnDictT& slnDict, countDictT& countDict, const size_t& k, const string& outputName ) {

        // Compute the *weighted* probability of each edge being in
        // the top k distinct scoring solutions

        // Each leaf has a single solution which is, by definition, of optimal cost
        for ( const auto& vit : order ) {
            if ( H->incident(vit).size() == 0 ) {countDict[vit].push_back( make_tuple(slnDict[vit][0].cost, cl_I(1)) );}
        }

        // For each edge, count the number of solutions having each score
        unordered_map< size_t, unordered_map< double, cl_I > > edgeCountMap;
        unordered_map< size_t, unordered_map< double, double > > edgeProbMap;

        // A map holding which edges are used to obtain *an* optimal
        // solution for each vertex
        Google<size_t, unordered_map<double, unordered_set<size_t> >>::Map usedEdges;
        usedEdges.set_empty_key(std::numeric_limits<size_t>::max());

        auto N = H->order();
        size_t ctr=0;

        // For each vertex, in topological order (backward)
        for ( auto vit = order.begin(); vit != order.end(); ++vit, ++ctr ) {
            cerr << "\r\rProcessed " << 100.0*(static_cast<float>(ctr)/N) << "% of the vertices";
            auto vert = H->vertex(*vit);

            // The cost for deriving this vertex that a new derivation has to beat
            double cval = ( slnDict.find(*vit) == slnDict.end() ) ? std::numeric_limits<double>::max() : slnDict[*vit][0].cost;

            map< double, vector<tuple<size_t,cl_I> > > edgeCostMap;
            // Soft constraint for edge existence
            auto ePen = existencePenalty(ti, vert, penalty);

            // loop over all incoming edges and compute the # of
            // solutions over each edge as well as that solution's cost
            for ( const auto& e : H->incident(*vit) ) {
                auto edge = H->edge(e);
                auto hind = edge.head();
                auto w = edge.weight();
                auto tvert = H->vertex(edge.tail()[0]);

                auto tval = ePen + w;
                cl_I numEdgeSln(1);

                auto currentEdgeSlns = countEdgeSolutions( tval, edge.tail(), countDict, k );

                for ( const auto& ent : currentEdgeSlns ) {
                    double score; cl_I count;
                    tie(score, count) = ent;
                    auto insertIt = edgeCostMap.find( score );
                    auto edgeContrib = make_tuple(e, count);

                    if ( insertIt == edgeCostMap.end() ) {
                        vector< tuple<size_t, cl_I > > r(1, edgeContrib);
                        edgeCostMap[ score ] = r;
                    } else {
                        insertIt->second.push_back( edgeContrib );
                    }
                    edgeCountMap[ e ][ score ] = count;
                }
            }
            // If we traversed any optimal edges
            if ( edgeCostMap.size() > 0 ) { //edgeCosts.size() > 0 ) {
                //cerr << "AM I CRAZY?\n";
                typedef tuple<double, cl_I, size_t> edgeSlnT;
                double minCost = std::numeric_limits<double>::max();
                size_t mk = std::min( edgeCostMap.size(), k );
                size_t ctr = 0;

                for ( auto cmIt = edgeCostMap.begin(); cmIt != edgeCostMap.end() && ctr < mk; ++cmIt, ++ctr ) {
                    auto score = cmIt->first;
                    //cerr << "ctr = " << ctr << ", score = " << score << "\n";
                    const auto& providingEdges = cmIt->second;
                    minCost = std::min(minCost, score);
                    cl_I numSln(0);
                    for ( const auto& edgeCount : providingEdges ) {
                        size_t edgeInd; cl_I count;
                        tie(edgeInd, count) = edgeCount;
                        if ( usedEdges[*vit].find(score) == usedEdges[*vit].end() ) {
                            usedEdges[*vit][score] = { edgeInd };
                        } else {
                            usedEdges[*vit][score].insert( edgeInd );
                        }
                        numSln += count;
                    }
                    countDict[*vit].push_back( make_tuple(score, numSln) );

                    for ( const auto& edgeCount : providingEdges ) {
                        size_t edgeInd; cl_I count;
                        tie(edgeInd, count) = edgeCount;
                        double edgeProb = double_approx( count / numSln );
                        edgeProbMap[edgeInd][score] = edgeProb;
                    }

                }



                // Find the minimum cost edge
                Google<Derivation::flipT>::Set fs;
                fs.set_empty_key( Derivation::flipT(-1,-1,"") );
                slnDict[*vit][0] = Derivation( minCost, 0, vector<size_t>(), fs);

            }
        }
        // loop over verts
        cerr << "\n";


        typedef Google< size_t, double >::Map probMapT;

        auto getOrElse = [] ( probMapT& pm, const size_t& key, double alt ) {
            return (pm.find(key) == pm.end()) ? alt : pm[key];
        };

        cerr << "Backward step \n";

        probMapT probMap;
        probMap.set_empty_key( std::numeric_limits<size_t>::max() );

        FlipKey rootKey( t->getRootId(), t->getRootId(), false, false);
        auto rootInd = H->index(rootKey);
        // The root gets a probability of 1
        probMap[rootInd] = 1.0;//cl_RA(1);
        //probMap[rootInd] = cl_RA(1);

        ctr = 0;
        size_t tot = order.size();

        // Compute the probabilities (forward)
        // Loop over all vertices in reverse topological order
        for ( auto vit = order.rbegin(); vit != order.rend(); ++vit, ++ctr ) {
            cerr << "\r\rprocessing node " << ctr << "/" << tot;

            auto key = H->vertex(*vit);
            auto parentProb = probMap[*vit];

            auto alphas = computeAlphasDouble( countDict[*vit] );

            // for each top-k score at this node
            for ( size_t i = 0; i < countDict[*vit].size(); ++i ) {
                double pScore; cl_I pCount;
                tie(pScore, pCount) = countDict[*vit][i];

                // for all incoming edges contributing to this score
                for ( const auto& e : usedEdges[*vit][pScore] ) {
                    auto condProb = edgeProbMap[e][pScore];//cl_RA(edgeCountMap[e][pScore] / pCount);
                    // for all tail vertices of this edge
                    auto tail = H->getTail(e);
                    for ( const auto& tind : tail ) {
                        //auto currentProb = getOrElse( probMap, tind, cl_RA(0) );
                        probMap[tind] += (parentProb * ( alphas[i] * condProb ));
                    }
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



    void viterbi( unique_ptr<ForwardHypergraph>& H, Tree* t, TreeInfo& ti, double penalty, const vector<size_t>& order, slnDictT& slnDict ) {
        auto N = H->order();
        size_t ctr=0;
        for ( auto vit = order.begin(); vit != order.end(); ++vit, ++ctr ) {
            cerr << "\r\rProcessed " << 100.0*(static_cast<float>(ctr)/N) << "% of the vertices";
            auto vert = H->vertex(*vit);
            //cout << "Processing : [" << t->getNodeName(vert.u()) << ", " << t->getNodeName(vert.v()) << " : (" << vert.f() << vert.r() << ")]\t";
            double cval = std::numeric_limits<double>::max();
            if ( slnDict.find(*vit) != slnDict.end() ) {
                cval = slnDict[*vit][0].cost;
            }

            //    cout << "COST = " << cval << "\n";

            auto existencePenalty = 0.0;
            if (vert.f() || vert.r()) {
                if (ti.intervalDistance( vert.u(), vert.v() ) > 0.0) {
                    existencePenalty = penalty * ti.intervalDistance( vert.u(), vert.v() );
                }
            }

            for ( auto& e : H->incident(*vit) ) {
                auto edge = H->edge(e);
                auto hind = edge.head();
                auto w = edge.weight();
                auto tvert = H->vertex(edge.tail()[0]);

                //if ( w > 0 ) { cout << "WEIGHt > 0 "; }
                auto tval = w;

                Google<Derivation::flipT>::Set flips;
                flips.set_empty_key( make_tuple(-1,-1,""));
                //set<Derivation::flipT> flips;
                string children;
                for ( auto tn : edge.tail() ){
                    auto cvert = H->vertex(tn);
                    tval += slnDict[tn][0].cost;
                    /*std::set_union( slnDict[tn][0].flips.begin(), slnDict[tn][0].flips.end(),
                                    slnDict[tn][0].flips.begin(), slnDict[tn][0].flips.end(),
                                    std::inserter( flips, flips.begin() ) );
                    */
                    for ( auto& f : slnDict[tn][0].flips ) { flips.insert(f); }

                    //children += "[ " + t->getNodeName(cvert.u()) + ", " + t->getNodeName(cvert.v()) + " : (" + lexical_cast<string>(cvert.f()) + lexical_cast<string>(cvert.r()) +")] ";
                }

                auto ft = flipType(vert,tvert);
                if ( ft != "n" ) { flips.insert( make_tuple(vert.u(), vert.v(), ft ) ); }

                if (tval <= cval) {
                    vector<size_t> bp( edge.tail().size(), 0 );
                    slnDict[hind][0] = Derivation(tval, e, bp, flips);
                    //cout << " updated cost from " << cval << " -> " << tval << " " << ((w > 0) ? "by flipping and" : "") <<  " using " << children << "\n";
                    cval = tval;
                }
            }
        }
        cerr << "\n";
    }




    /** Alg 3 from the paper */
    typedef std::pair<size_t, Derivation> TaggedDerivT;

    Google< tuple<size_t,size_t> >::Set recursedStore;
    Google<size_t, vector<Derivation>*>::Map cand;
    //unordered_map< size_t, vector<Derivation>* > cand;

    void initKBest() {
        recursedStore.set_empty_key( make_tuple( std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max() ) );
        cand.set_empty_key( std::numeric_limits<size_t>::max() );
    }

    template<typename T>
    void printVector( const vector<T>& v ) {
        cerr << "[";
        for ( auto& e : v ) { cerr << e << " "; } cerr << "]";
    }

    void getCandidates( const unique_ptr<ForwardHypergraph>& H, TreeInfo& ti, double penalty, slnDictT& d, size_t v, size_t k, vector<Derivation>* kbest ){
        // There is only one candidate for any leaf node
        auto incoming = H->incident(v);
        auto numEdges = incoming.size();
        if ( numEdges == 0 ) { return; }
        // For internal nodes
        // List of potential candidates
        vector<Derivation> temp;
        // The head vertex and it's constituent tree nodes
        auto hvert = H->vertex(v);
        auto a = hvert.u(); auto b = hvert.v();

        auto existencePenalty = 0.0;
        if (hvert.f() || hvert.r()) {
            if (ti.intervalDistance( a, b ) > 0.0) {
                existencePenalty = penalty * ti.intervalDistance( a, b );
            }
        }


        // Initalize the heap with the 1-best solution for each incoming edge
        for ( auto e : incoming ) {
            // Compute the score of this derivation
            auto edge = H->edge(e);
            auto w = existencePenalty + edge.weight();
            auto tail = edge.tail();
            // The DBP pointing to the 1-best solution of all of tail vertices of this edge
            vector<size_t> j( tail.size(), 0);
            // A Tail vertex for this edge
            auto tvert = H->vertex( tail[0] );
            auto ft = flipType(hvert, tvert);

            //cerr << "candid deriv : "; printVector( j ); cerr << "\n";

            // The weight of this derivation
            auto s = w;
            // The set of effective flips for this derivation
            Google<Derivation::flipT>::Set ed;
            ed.set_empty_key( make_tuple(-1,-1,""));
            // set<Derivation::flipT> ed;

            for ( auto t : tail ){
                s += d[t][0].cost;
                /*std::set_union( d[t][0].flips.begin(), d[t][0].flips.end(),
                                d[t][0].flips.begin(), d[t][0].flips.end(),
                                std::inserter( ed, ed.begin() ) );
                */
                for ( auto& f : d[t][0].flips ) { ed.insert(f); }
            }

            if ( ft != "n" ) { ed.insert( make_tuple(a,b,ft) ); }

            //make_heap( temp.begin(), temp.end() );
            temp.push_back( Derivation(s, e, j, ed) );
            //push_heap( temp.begin(), temp.end() );

        }
        // We only need the top k solutions here
        // kbest = heapq.nsmallest( max(k,numEdges), temp );
        // Heapify the list and return it
        k = std::min( temp.size(), std::max( k, numEdges ) );
        partial_sort( temp.begin(), temp.begin()+k, temp.end(),
                      [] ( const Derivation& d1, const Derivation& d2 ) { return d1.cost < d2.cost; } );
        kbest->reserve(k);
        copy( temp.begin(), temp.begin()+k, std::back_inserter(*kbest) );
        //for ( auto c = temp.begin(); c != temp.begin()+k; ++c ) {
        //    kbest.push_back(*c);
        //}
        //cerr << "|temp| = " << temp.size() << ", |kbest| = " << kbest.size() << "\n";
    }

    bool sameEffect( const Derivation& d0, const Derivation& d1 ) {
        if ( d0.cost != d1.cost ) { return false; }
        // otherwise, same scores, check the flips
        // vector<Derivation::flipT> res;
        // std::set_symmetric_difference( d0.flips.begin(), d0.flips.end(), d1.flips.begin(), d1.flips.end(), std::back_inserter(res) );
        // return (res.size() == 0) ? true : false;

        // for ( auto f0 : d0.flips ) {
        //   if ( d1.flips.find(f0) == d1.flips.end() ) { return false; }
        // }

        bool hasDiffElement = any_of( d0.flips.begin(),
                                      d0.flips.end(),
                                      [&] ( const tuple<int,int,string>& d ) { return d1.flips.find(d) == d1.flips.end(); } );

        // for ( auto f1 : d1.flips ) {
        //     if ( d0.flips.find(f1) == d0.flips.end() ) { return false; }
        // }
        //cerr << "\n\nSame Effect:\n" << d0 << "\n" << d1 << "\n\n";
        return !hasDiffElement;
    }

    void lazyKthBest( unique_ptr<ForwardHypergraph>& H, TreeInfo& ti, double penalty, slnDictT& d, size_t v, size_t k, size_t kp );

    void lazyNext( unique_ptr<ForwardHypergraph>& H, TreeInfo& ti, double penalty, slnDictT& d, vector<Derivation>* locCand, size_t eind,
               const vector<size_t>& j, size_t kp) {

        auto edge = H->edge(eind);
        auto tail = edge.tail();
        auto hvert = H->vertex(edge.head());
        auto tvert = H->vertex(tail[0]);
        auto a = hvert.u(); auto b = hvert.v();

        auto ft = flipType(hvert, tvert);

        auto existencePenalty = 0.0;
        if (hvert.f() || hvert.r()) {
            if (ti.intervalDistance( a, b ) > 0.0) {
                existencePenalty = penalty * ti.intervalDistance( a, b );
            }
        }

        // For each tail node in this edge
        size_t i = 0;
        for ( auto vind : tail ) {
            //cerr << "before : deriv ";
            //printVector(j);
            //cerr << "\n";

            // Ask for the next best derivation at this tail node
            vector<size_t> jp(j); jp[i] += 1;

            //cerr << "after : deriv ";
            //printVector(jp);
            //cerr << "\n";

            // Compute it recursively
            //cerr << "Calling lazyKthBest on " << vind << "  ";
            //printVector(jp);
            lazyKthBest( H, ti, penalty, d, vind, jp[i], kp );



            // If this derivation exists
            if (jp[i] < d[vind].size() ) {
                // Cost of the derivation -- edge cost + sum of children costs
                // Effective edges -- Union of the effective edges of
                // the children

                auto w = edge.weight();
                auto s = existencePenalty + w;
                Google<Derivation::flipT>::Set ed;
                ed.set_empty_key( make_tuple(-1,-1,""));
                //set<Derivation::flipT> ed;
                /*
                for (size_t ind = 0; ind < jp.size(); ind+=2 ) {
                    // Process a pair
                    if ( jp.size() - ind > 1 ) {
                        auto vp1 = tail[ind]; auto k1 = jp[ind];
                        auto vp2 = tail[ind+1]; auto k2 = jp[ind+1];

                        s += (d[vp1][ k1].cost + d[vp2][k2].cost);

                        std::set_union( d[vp1][k1].flips.begin(), d[vp1][k1].flips.end(),
                                        d[vp2][k2].flips.begin(), d[vp2][k2].flips.end(),
                                        std::inserter( ed, ed.begin() ) );

                    } else { // Process a single element
                        auto vp1 = tail[ind]; auto k1 = jp[ind];
                        s += d[vp1][ k1].cost;
                        std::set_union( d[vp1][k1].flips.begin(), d[vp1][k1].flips.end(),
                                        d[vp1][k1].flips.begin(), d[vp1][k1].flips.end(),
                                        std::inserter( ed, ed.begin() ) );


                    }
                }*/


                for (size_t ind = 0; ind < jp.size(); ++ind ) {
                    auto vp = tail[ind]; auto k = jp[ind];
                    s += d[ vp ][ k ].cost;
                    for ( auto& f : d[vp][k].flips ) { ed.insert(f); }
                }


                if ( ft != "n" ) {
                    // If this edge is effective, then it is added to the effective edges of this derivation
                    ed.insert( make_tuple(a,b,ft) );
                }
                // If this derivation is not in the heap, then add
                // it

                Derivation deriv(s, eind, jp, ed);
                auto foundDer = false;
                for ( auto& od : *locCand ) {
                    if ( sameEffect(deriv, od) ) { foundDer = true; break; }
                }
                if ( !foundDer ) {
                    locCand->push_back( deriv );
                    push_heap( locCand->begin(), locCand->end(),
                               [] ( const Derivation& d0, const Derivation& d1) { return d0.cost > d1.cost; } );
                }

            }

            i++;
        }
    }


    void lazyKthBest( unique_ptr<ForwardHypergraph>& H, TreeInfo& ti, double penalty, slnDictT& d, size_t v, size_t k, size_t kp ) {
        // If this is our first vist to vertex v, then
        // populate its candidate list

        //cerr << "Searching for candidates for " << v << "\n";
        if ( cand.find(v) == cand.end() ) {
            //cerr << "There were none; adding to the cand map \n";
            //if ( kp == 0 ) { recursedStore.set_empty_key( make_tuple(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()) ); }
            cand[v] = new vector< Derivation >();
            //cerr << "done \n";
            getCandidates(H, ti, penalty, d, v, kp, cand[v]);
            /*
            cerr << "Added candidates for node " << v << "\n[\n";
            for ( auto& c : cand[v] ) {
                cerr << c << "\n";
            } //cerr << "]\n";
            */
        }

        // Until we have the required number of derivations
        while ( d[v].size() <= k ) {
            // If we've already computed the last derivation
            auto lastInd = 0;
            if ( d[v].size() > 0 ) {
                // Get this derivation
                auto lastInd = max_element( d[v].begin(), d[v].end(),
                                            [] ( const TaggedDerivT& a, const TaggedDerivT& b ) { return a.first < b.first; } )->first;
                auto deriv = d[v][lastInd];
                //auto deriv = d[v].back();

                // Update the heap adding the last derivation's
                // successor
                if ( recursedStore.find( make_tuple(v,lastInd) ) == recursedStore.end() ) {
                    //cerr << "Calling LazyNext from branch 1\n";
                    //cerr << "lastInd = " << lastInd << " deriv = " << deriv << "\n";
                    if ( H->incident(v).size() > 0 ) {
                        //cerr << " current bp = ";
                        //printVector(bp);
                        //cerr << "calling lazy next on cand[" << v << "]\n";
                        lazyNext(H, ti, penalty, d, cand[v], deriv.target, deriv.bp, kp);
                    }
                    //if (v, lastInd) not in recursedStore: recursedStore[ (v, lastInd) ] = False
                    recursedStore.insert( make_tuple(v, lastInd) );
                } else {
                    // Do Nothing
                    // cerr << "I've already looked for the successors of " << d[v][lastInd]  << "\n";
                    //cerr << "There are " << cand[v].size() << " candidates\n";
                }
            }
            // Get the next best candidate from the heap
            if ( cand[v]->size() > 0 ) {
                //cerr << "Examining candidates " << cand[v].size() << "\n";
                //cerr << "popping candidate from heap\n";
                pop_heap( cand[v]->begin(), cand[v]->end(),
                          [] ( const Derivation& d0, const Derivation& d1) { return d0.cost > d1.cost; } );
                auto deriv = cand[v]->back();
                //cerr << "done; deriv is " << deriv << " \n";
                //cerr << "before pop back, cand[ " << v << "].size() is " << cand[v]->size() << "\n";
                cand[v]->pop_back();
                //cerr << "really done\n";

                //cerr << "candidate is : " << deriv << "\n";
                /*
                auto lastInd = (max_element( d[v].begin(), d[v].end(),
                                             [] ( const TaggedDerivT& a, const TaggedDerivT& b ) { return a.first < b.first; } ))->first;
                */
                auto existingElem = find_if( d[v].begin(), d[v].end(), [&] ( const TaggedDerivT& od ) { return sameEffect(deriv,od.second); } );
                auto found = (existingElem != d[v].end());

                // Always add a derivation if none exist
                if ( d[v].size() == 0 || ! found ) {
                    d[v][ d[v].size() ] = deriv;
                    //d[v].push_back(deriv);

                } else {
                    // If there is no derivation with this same score and set of effective flips, then add it
                    //cerr << "skipping solution " << deriv << "\n";
                    if ( H->incident(v).size() > 0 ) {
                        //cerr << "2 calling lazyNext on cand[" << v << "]\n";
                        lazyNext(H, ti, penalty, d, cand[v], deriv.target, deriv.bp, kp);
                    }
                }
            } else {
                // We have no more derivations for this node!
                //cerr << "there are 0 candidates for " << H->vertex(v) << "\n";
                return;
            }
        }
    }


}

#endif // MULTI_OPT_HPP
