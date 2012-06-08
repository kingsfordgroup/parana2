#include "MultiOpt.hpp"
#include "ProgressDisplay.hpp"
#include <cmath>
#include <utility>
#include <boost/timer/timer.hpp>
#include <boost/heap/fibonacci_heap.hpp>


/** Google's dense hash set and hash map **/
#include <google/dense_hash_set>
#include <google/dense_hash_map>
#include <google/sparse_hash_map>


namespace MultiOpt {

    const auto none = make_tuple(false,false);
    const auto fw = make_tuple(true,false);
    const auto rev = make_tuple(false,true);
    const auto both = make_tuple(true,true);

    static unordered_map<tuple<bool,bool>, string> flipStrMap = {
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
        auto none = make_tuple(false,false);
        auto fw = make_tuple(true,false);
        auto rev = make_tuple(false,true);
        auto both = make_tuple(true,true);

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
        auto none = make_tuple(false,false);
        auto fw = make_tuple(true,false);
        auto rev = make_tuple(false,true);
        auto both = make_tuple(true,true);

        costMapFunT costDict;
        costDict[ none ][ none ] = [=](const double& pf, const double& pr) {
            return make_tuple( 0.0, "n" );
        };
        costDict[ none ][ fw ] = [=](const double& pf, const double& pr) {
            return make_tuple( pf*cc, "f+" );
        };
        costDict[ none ][ rev ] = [=](const double& pf, const double& pr) {
            return make_tuple( pr*cc, "r+" );
        };
        costDict[ none ][ both ] = [=](const double& pf, const double& pr) {
            return make_tuple( undirected ? pf*cc : (pf+pr)*cc, "b+" );
        };


        costDict[ rev ][ none ] = [=](const double& pf, const double& pr) {
            return make_tuple( pr*dc, "r-" );
        };
        costDict[ rev ][ fw ] = [=](const double& pf, const double& pr) {
            return make_tuple( pf*cc + pr*dc, "f+r-" );
        };
        costDict[ rev ][ rev ] = [=](const double& pf, const double& pr) {
            return make_tuple( (1.0-pr)*dc, "n");
        };
        costDict[ rev ][ both ] = [=](const double& pf, const double& pr) {
            return make_tuple( pf*cc + (1.0-pr)*dc , "f+");
        };

        costDict[ fw ][ none ] = [=](const double& pf, const double& pr) {
            return make_tuple( pf*dc, "f-");
        };
        costDict[ fw ][ fw ] = [=](const double& pf, const double& pr) {
            return make_tuple( (1.0-pf)*dc, "n");
        };
        costDict[ fw ][ rev ] = [=](const double& pf, const double& pr) {
            return make_tuple( pf*dc + pr*cc, "f-r+");
        };
        costDict[ fw ][ both ] = [=](const double& pf, const double& pr) {
            return make_tuple( (1.0-pf)*dc + pr*cc, "r+");
        };

        costDict[ both ][ none ] = [=](const double& pf, const double& pr) {
            return make_tuple( undirected ? dc : (pf+pr)*dc, "b-" );
        }; // was pf*dc
        costDict[ both ][ fw ] = [=](const double& pf, const double& pr) {
            return make_tuple( (1.0-pf)*dc + pr*dc, "r-" );
        };
        costDict[ both ][ rev ] = [=](const double& pf, const double& pr) {
            return make_tuple( pf*dc + (1.0-pr)*dc, "f-" );
        };
        costDict[ both ][ both ] = [=](const double& pf, const double& pr) {
            return make_tuple( undirected ? (1.0-pf)*dc : (1.0-pf)*dc + (1.0-pr)*dc, "n" );
        };

        return costDict;
    }

    selfCostMapT getSelfLoopCostDict ( double cc, double dc, bool directed ) {
        auto undirected = ! directed;
        auto none = make_tuple(false,false);
        auto fw = make_tuple(true,false);
        auto rev = make_tuple(false,true);
        auto both = make_tuple(true,true);

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
        auto none = make_tuple(false,false);
        auto fw = make_tuple(true,false);
        auto rev = make_tuple(false,true);
        auto both = make_tuple(true,true);

        selfCostMapFunT selfLoopCostDict;
        selfLoopCostDict[ none ][ true ] = [=] (const double& p) {
            return make_tuple( p*cc, "b+" );
        };
        selfLoopCostDict[ none ][ false ] = [=] (const double& p) {
            return make_tuple( p*dc, "n" );
        };

        selfLoopCostDict[ rev ][ true ] = [=] (const double& p) {
            return make_tuple( (1.0-p)*dc, "n" );
        };
        selfLoopCostDict[ rev ][ false ] = [=] (const double& p) {
            return make_tuple( p*dc, "b-" );
        };

        selfLoopCostDict[ fw ][ true ] = [=] (const double& p) {
            return make_tuple( (1.0-p)*dc, "n" );
        };
        selfLoopCostDict[ fw ][ false ] = [=] (const double& p) {
            return make_tuple( p*dc, "b-" );
        };

        selfLoopCostDict[ both ][ true ] = [=] (const double& p) {
            return make_tuple( (1.0-p)*dc, "n" );
        };
        selfLoopCostDict[ both ][ false ] = [=] (const double& p) {
            return make_tuple( dc, "b-" );
        }; // was p*dc
        return selfLoopCostDict;
    }

    string flipType( const FlipKey& hvert, const FlipKey& tvert ) {
        return flipDict[ make_tuple(hvert.f(), hvert.r())][ make_tuple(tvert.f(), tvert.r()) ];
    }

    template<typename GT>
    std::unordered_set<int> projectToReversedGraph( unique_ptr<ForwardHypergraph>& H, GT& G ) {
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

    void topologicalOrder( unique_ptr<ForwardHypergraph>& H, vector<size_t>& order ) {
        using boost::adjacency_list;
        using boost::vecS;
        using boost::directedS;

        typedef adjacency_list<vecS, vecS, directedS> graphT;
        graphT G( H->order() );

        vector<size_t> torder;
        std::unordered_set<int> vset = projectToReversedGraph( H, G );
        boost::topological_sort(G, std::back_inserter(torder));

        // Remove all of the vertices from "order" if they are not in vset
        cerr << "ORDER BEFORE = " << torder.size() << "\n";
        std::copy_if(torder.begin(), torder.end(), std::back_inserter(order), [=](const int& id) -> bool { return vset.find(id) != vset.end(); } );
        //std::copy()
        cerr << "ORDER AFTER = " << order.size() << "\n";
    }

    /**
     *  Compute the penalty for this edge to exist based on difference
     *  between the existence intervals of the endpoints and the
     *  penalty factor.
     */
    template <typename T>
    double existencePenalty( const TreeInfo& ti, const T& vert, const double& penalty, const double& travWeight ) {
        if ( travWeight > 0 ) {
            auto dist = ti.intervalDistance( vert.u(), vert.v() );
            if ( dist > 0.0 ) {
                return penalty * dist;
            }
        }
        return 0.0;
    }

    unique_ptr<ForwardHypergraph>  buildSolutionSpaceGraph( const TreePtrT& t,
            const TreeInfo& ti,
            double cc,
            double dc,
            double penalty,
            bool directed ) {
        boost::timer::auto_cpu_timer timer;

        Google<int>::Set leafSet;
        leafSet.set_empty_key(-1);
        for ( auto l : t->getLeavesId() ) {
            leafSet.insert(l);
        }
        auto isLeaf = [&]( const int& nid ) -> bool { return leafSet.find(nid) != leafSet.end(); };//t->isLeaf(nid); };
        auto isInternal = [&]( const int& nid ) -> bool { return leafSet.find(nid) == leafSet.end(); };
        auto isLost = [&]( int nid ) -> bool { return (t->getNodeName(nid)).find("LOST") != std::string::npos; };

        auto rootId = t->getRootId();
        auto rootName = t->getNodeName( rootId );

        auto fauxRoot = "#preroot#";

        unique_ptr<ForwardHypergraph> slnSpaceGraph( new ForwardHypergraph() );

        costMapT costMap( getCostDict(cc,dc,directed) );
        selfCostMapT selfLoopCostMap( getSelfLoopCostDict(cc,dc,directed) );

        auto addIncomingHyperedge = [=,&costMap,&selfLoopCostMap,&slnSpaceGraph,&penalty,&t] (const FlipKey& k, const int& rnode, const int& onode ) {
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
                        auto w = get<0>(selfLoopCostMap[ make_tuple(k.f(), k.r()) ][ dualFlipLL.f() || dualFlipLL.r() ]);
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
                        auto w = get<0>(costMap[ make_tuple(k.f(),k.r()) ][ make_tuple(dualFlipL.f(), dualFlipL.r()) ]);//directed ? 2.0 : 1.0;//2.0 ? directed : 1.0;
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
                            auto w = get<0>(costMap[ make_tuple(k.f(), k.r() ) ][ make_tuple(fwFlipL.f(), fwFlipL.r()) ]);
                            w += existencePenalty(ti, k, penalty, w);
                            if ( std::isfinite(w) ) {
                                slnSpaceGraph->addEdge( fwFlip, k , w );
                            }
                        }
                        if ( revFlip.size() > 0 ) {
                            auto w = get<0>(costMap[ make_tuple(k.f(), k.r() ) ][ make_tuple(revFlipL.f(), revFlipL.r()) ]);
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
                        (!differentExtantNetworks(ti,u,v) &&
                         !(ti.inSubnodesOf(u,v) ||  ti.inSubnodesOf(v,u)) &&
                         !(isLost(u) || isLost(v)) )) {
                    ctr1 += 1;
                    slnSpaceGraph->addVertex( FlipKey( u, v, false, false ) );
                    if ( ! t->isRoot(v) ) {
                        ctr2 += 1;
                        slnSpaceGraph->addVertex( FlipKey( u, v, true, true ) );
                    }
                    if ( directed && u!=v ) {
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

        auto vstr = [&]( const FlipKey& vert ) -> string {
            auto uname = t->getNodeName(vert.u());
            auto vname = t->getNodeName(vert.v());
            if (uname > vname) {
                auto tmp = vname;
                vname = uname;
                uname = tmp;
            }
            auto fstr = vert.f() ? "true" : "false";
            auto rstr = vert.r() ? "true" : "false";
            return uname+"\t"+vname+"\t"+fstr+"\t"+rstr;
        };

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
        cerr << "Hypergraph size = " << slnSpaceGraph->size() << "\n";
        return slnSpaceGraph;
    }


    template< typename GT >
    void leafCostDict( unique_ptr<ForwardHypergraph>& H, TreePtrT& T, GT& G, bool directed, double cc, double dc, slnDictT& slnDict ) {
        /*
          Given the duplication tree T, the root vertex rv, the extant graph G and
          the constraints, fill in the values for the leaf nodes of the hypergraph
        */
        typedef typename boost::graph_traits< GT >::edge_descriptor EdgeT;
        boost::timer::auto_cpu_timer timer;
        auto undirected = !directed;
        auto isLost = [&]( int nid ) -> bool { return (T->getNodeName(nid)).find("LOST") != std::string::npos; };
        // Cost of going from the state encoded by a node to the state of a true graph edge

        auto none = make_tuple(false,false);
        auto fw = make_tuple(true,false);
        auto rev = make_tuple(false,true);
        auto both = make_tuple(true,true);

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
        for( size_t i = 0; i < N; ++i) {
            auto elist = H->incident(i);
            if ( elist.size() == 0 ) {
                leafHypernodes.push_back(i);
            }
        }

        //cerr << "# of leaf hypernodes = " << leafHypernodes.size() << "\n";

        typedef typename GT::vertex_descriptor NodeT;
        typedef unordered_set<NodeT> NodeSetT;

        // Is the node e contained in the set s?
        auto contains = [] ( const NodeSetT& s, NodeT e ) {
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
        for( auto n : leafHypernodes ) {
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
                es.set_empty_key( make_tuple(-1,-1,""));
                nlost += 1;
                slnDict[n] = { {0,Derivation(lostCost, n, ev, es)} };
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
                    tie(fedge, d_f) = edge(u,v,G);
                    double w_f = d_f ? G[fedge].weight : 0.0;
                    tie(redge, d_r) = edge(v,u,G);
                    double w_r = d_r ? G[redge].weight : 0.0;
                    if ( undirected ) {
                        assert( w_f == w_r );
                    }

                    tweight += (w_f + w_r) * (d_f + d_r) * (nd.f() + nd.r());

                    //auto costFlipProb = costDict[ make_tuple(f,r) ][ make_tuple(d_f,d_r) ];
                    auto costFlipProb = costFunDict[ make_tuple(f,r) ][ make_tuple(d_f,d_r) ](w_f, w_r);

                    /*
                    if (costFlip != costFlipProb) {
                        cerr << "whoops for transition (" << f << ", " << r << ") => (" << d_f << ", " << d_r << "), and (w_f, w_r) = (" << w_f << ", " << w_r << ")\n";
                        cerr << "costFlip = (" << get<0>(costFlip) << ", " << get<1>(costFlip) << "), but costFlipProb = (" << get<0>(costFlipProb) << ", " << get<1>(costFlipProb) << ")\n";
                        exit(1);
                    }
                    */

                    auto cost = get<0>(costFlipProb);
                    auto flip = get<1>(costFlipProb);

                    Google<Derivation::flipT>::Set effectiveEdges;
                    effectiveEdges.set_empty_key( make_tuple(-1,-1,""));

                    if ( flip != "n" ) {
                        nef += 1;
                        effectiveEdges.insert( make_tuple(nd.u(),nd.v(),flip) );
                    }
                    vector<size_t> ev;
                    slnDict[n] = { {0,Derivation(cost, n, ev, effectiveEdges)} };

                } else {
                    EdgeT e;
                    bool hasSelfLoop;
                    tie(e, hasSelfLoop) = edge(u,v,G);
                    double w_l = hasSelfLoop ? G[e].weight : 0.0;

                    tweight += w_l * (hasSelfLoop) * (nd.f() + nd.r());
                    //auto costFlipProb = selfLoopCostDict[ make_tuple(f,r) ][ hasSelfLoop ];
                    auto costFlipProb = selfLoopCostFunDict[ make_tuple(f,r) ][ hasSelfLoop ]( w_l );
                    /*
                    if (costFlip != costFlipProb) {
                        cerr << "whoops for self loop transition (" << f << ", " << r << ") => (" << hasSelfLoop << "), and (w_l) = (" << w_l << ")\n";
                        cerr << "costFlip = (" << get<0>(costFlip) << ", " << get<1>(costFlip) << "), but costFlipProb = (" << get<0>(costFlipProb) << ", " << get<1>(costFlipProb) << ")\n";
                        exit(1);
                    }
                    */

                    auto cost = get<0>(costFlipProb);
                    auto flip = get<1>(costFlipProb);

                    Google<Derivation::flipT>::Set effectiveEdges;
                    effectiveEdges.set_empty_key( make_tuple(-1,-1,""));

                    if ( flip != "n" ) {
                        nef +=1;
                        effectiveEdges.insert( make_tuple(nd.u(),nd.v(),flip) );
                    }
                    vector<size_t> ev;
                    slnDict[n] = { {0,Derivation(cost, n, ev, effectiveEdges)} };

                } // ( u != v )
            } // ( lostU || lostV )
        } // loop over leaf hypernodes

        double sum = 0.0;
        //cerr << "slnDict size = " << slnDict.size() << "\n";
        for( auto n : slnDict ) {
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

    tuple<double, cl_I> getCostCount( unordered_map<size_t, vector<CostClass> >& tkd,
                                      const vector<size_t>& bp,
                                      const size_t& eid,
                                      unique_ptr<ForwardHypergraph>& H ) {
        auto edge = H->edge(eid);
        double cost = edge.weight();
        cl_I count(1);
        size_t i = 0;
        for( auto& tailNode : edge.tail() ) {
            // the cost class index
            size_t cci = bp[i];
            // We sum the costs and multiply the counts
            cost += tkd[tailNode][cci].cost();
            count *= tkd[tailNode][cci].total();
            i++;
        }
        return make_tuple(cost, count);
    }


    double getCost( unordered_map<size_t, vector<CostClass> >& tkd,
                    const vector<size_t>& bp,
                    const size_t& eid,
                    unique_ptr<ForwardHypergraph>& H ) {
        auto edge = H->edge(eid);
        double cost = edge.weight();
        size_t i = 0;
        for( auto& tailNode : edge.tail() ) {
            // the cost class index
            size_t cci = bp[i];
            // We sum the costs and multiply the counts
            cost += tkd[tailNode][cci].cost();
            i++;
        }
        return cost;
    }

    cl_I getCount( unordered_map<size_t, vector<CostClass> >& tkd,
                   const vector<size_t>& bp,
                   const size_t& eid,
                   unique_ptr<ForwardHypergraph>& H ) {
        auto edge = H->edge(eid);
        cl_I count(1);
        size_t i = 0;
        for( auto& tailNode : edge.tail() ) {
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
        const double& ecost,
        const vector<size_t>& tailNodes,
        countDictT& countDict,
        const size_t& k,
        bool printMe,
        unique_ptr<ForwardHypergraph>& H,
        TreePtrT& t ) {

        auto vstr = [&]( const FlipKey& vert ) -> string {
            auto uname = t->getNodeName(vert.u());
            auto vname = t->getNodeName(vert.v());
            if (uname > vname) {
                auto tmp = vname;
                vname = uname;
                uname = tmp;
            }
            auto fstr = vert.f() ? "true" : "false";
            auto rstr = vert.r() ? "true" : "false";
            return uname+"\t"+vname+"\t"+fstr+"\t"+rstr;
        };

        // product pointers
        std::vector< size_t > elemSizes;
        elemSizes.reserve(tailNodes.size());
        double cost = ecost;
        for ( const auto& t : tailNodes ) {
            elemSizes.push_back( countDict[t].size() );
            cost += get<0>(countDict[t].front());
        }

        vector<dvsT> pq(1, make_tuple(cost, vector<size_t>(tailNodes.size(), 0)));
        QueueCmp<dvsT> ord;


        std::function< double( const vector<size_t>& ) > computeScore = [&] ( const vector<size_t>& inds ) -> double {
            size_t numNodes = tailNodes.size();
            double cost = ecost;
            for ( size_t i = 0; i < numNodes; ++i ) {
                cost += get<0>(countDict[ tailNodes[i] ][ inds[i] ]);
            }
            return cost;
        };

        auto round3 = [=]( double num ) -> double {
            double result = num * 1000;
            result = std::floor(result);
            result = result / 1000;
            return result;
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
            std::tie(cost,inds) = pq.front();
            std::pop_heap( pq.begin(), pq.end(), ord );
            pq.pop_back();

            // Compute the number of ways we can obtain this solution
            cl_I numSlns(1);
            for ( size_t i = 0; i < inds.size(); ++i ) {
                if(printMe) {
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
                    for( auto es: edgeSlns)  {
                        cerr << get<1>(es) << ", ";
                    }
                    cerr << "\t new guy (" << cost << ", " << numSlns << ")\n";
                }
                edgeSlns.push_back( make_tuple( cost, numSlns ) );
            } else { // we found a solution of this score
                if (printMe) {
                    cerr << "IN ELSE: edgeSlns = ";
                    for( auto es: edgeSlns)  {
                        cerr << get<1>(es) << ", ";
                    }
                    cerr << "\t new guy (" << cost << ", " << numSlns << ")\n";
                }
                if( cost < get<0>(edgeSlns.back()) ) {
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
        if ( edgeSlns.size() == k+1 ) {
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

    vector<double> computeAlphasDouble( const vector<tuple<double, cl_I>>& slnVec, size_t k, const cl_I& total ) {
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
for (const auto& e : slnVec) {
                cerr << "score = " << get<0>(e) << ", count = " << get<1>(e) << "\n";
            }
            std::abort();
        }
        double diff = worstScore == bestScore ? 1.0 : worstScore - bestScore; // std::max(0.01, worstScore - bestScore);
        size_t N = slnVec.size();

        std::stringstream ss;
        ss << get<1>(slnVec.front());
        mpreal I = ss.str();
        ss.clear();
        ss << get<1>(slnVec.back());
        mpreal J = ss.str();
        ss.clear();

        //double scale = (mpfr::log( J ) - mpfr::log( I )).toDouble() / diff;
        //double scale = 2.0 * estimateThermodynamicBeta( slnVec, bestScore ); // (6.9*k) / diff; // (1.25 * k) / diff;
        //double scale = 1.8;
        double scale = 60.0 / diff;
        //double scale = (0.5 * k) / diff;//(2.0 * N) / diff;//(1.5*k) / diff;
        //std::cerr << " **** beta = " << scale << " **** \n";
        double sum(0.0);

for (const auto& e : slnVec) {
            double a = std::abs( (bestScore - get<0>(e)) * scale);
            scores.push_back( std::exp( -(a) ) );
            sum += scores.back();
        }

        double invSum = 1.0 / sum;
        vector<double> alphas;
        alphas.reserve(slnVec.size());
for (const auto& s : scores) {
            alphas.push_back( s * invSum );
        }
        return alphas;
    }

    void viterbiCount( unique_ptr<ForwardHypergraph>& H, TreePtrT& t, TreeInfo& ti, double penalty, const vector<size_t>& order,
                       slnDictT& slnDict, countDictT& countDict, const size_t& k,
                       const string& outputName, const vector<FlipKey>& outputKeys ) {

        // Compute the *weighted* probability of each edge being in
        // the top k distinct scoring solutions

        auto vstr = [&]( const FlipKey& vert ) -> string {
            auto uname = t->getNodeName(vert.u());
            auto vname = t->getNodeName(vert.v());
            if (uname > vname) {
                auto tmp = vname;
                vname = uname;
                uname = tmp;
            }
            auto fstr = vert.f() ? "true" : "false";
            auto rstr = vert.r() ? "true" : "false";
            return uname+"\t"+vname+"\t"+fstr+"\t"+rstr;
        };

        double costsum = 0.0;
        // Each leaf has a single solution which is, by definition, of optimal cost
        for ( const auto& vit : order ) {
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
for ( auto tup: ele.second ) {
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
        size_t ctr=0;

        auto oname = "EDGES.txt";
        std::fstream out( oname, std::fstream::out | std::fstream::trunc );

        // For each vertex, in topological order (up)
        for ( auto vit = order.begin(); vit != order.end(); ++vit, ++ctr ) {
            if ( !(ctr % 1000) || ctr == N-1 ) {
                cerr << "\r\rProcessed " << 100.0*(static_cast<float>(ctr)/N) << "% of the vertices";
            }
            auto vert = H->vertex(*vit);

            // Put the results in an ordered map -- the sorted
            // property will be useful later for maintaining only the
            // k-best score classes
            map< double, vector<tuple<size_t,cl_I> > > edgeCostMap;

            //cerr << "SLN FOR NODE " << vstr(vert) << "\n";
            out << vstr(vert);// << H->incident(*vit).size() << "\n";
            std::vector<cl_I> nedgesln;
            // loop over all incoming edges and compute the # of
            // solutions over each edge as well as that solution's cost
            for ( const auto& e : H->incident(*vit) ) {

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
                    for( auto el : currentEdgeSlns ) {
                        cerr << get<1>(el) << ", ";
                    }
                    cerr << "\n";
                }
                //if ( t->getNodeName(vert.u()) == "10866" )

                cl_I edgeSum(0);
                for ( const auto& ent : currentEdgeSlns ) {
                    edgeSum += get<1>(ent);
                }
                nedgesln.push_back(edgeSum);

                out << w << "\t" << currentEdgeSlns.size();
                size_t j = 0;
                for ( const auto& ent : currentEdgeSlns ) {
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
                for ( const auto& ent : currentEdgeSlns ) {
                    double score;
                    cl_I count;
                    tie(score, count) = ent;
                    auto edgeContrib = make_tuple(e, count);
                    edgeCostMap[score].push_back( edgeContrib );
                    edgeCountMap[ e ][ score ] = count;
                }
            }

            std::sort(nedgesln.begin(), nedgesln.end(),
            []( const cl_I& x, const cl_I& y) {
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
                    const auto& providingEdges = cmIt->second;
                    // the minimum cost incoming score
                    minCost = std::min(minCost, score);
                    // will count the number of solutions in this
                    // score class
                    cl_I numSln(0);

                    // Update the information at the derived vertices
                    for ( const auto& edgeCount : providingEdges ) {
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
                    for ( const auto& edgeCount : providingEdges ) {
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
                fs.set_empty_key( Derivation::flipT(-1,-1,"") );
                slnDict[*vit][0] = Derivation( minCost, 0, vector<size_t>(), fs);

            }
        }
        // loop over verts
        cerr << "\n";

        auto epname = "EPROBS.txt";
        std::fstream epout( epname, std::fstream::out | std::fstream::trunc );
        for( auto& em : edgeProbMap ) {
            //edgeProbMap.foreach{ case (eind, smap) =>
            auto eind = em.first;
            auto smap = em.second;
            epout << vstr( H->vertex(H->getHead(eind)) )  << "\t";
            epout << H->edge(eind).weight() << "\t" << smap.size()  << "\n";
            for( auto sp: smap ) {
                epout << sp.first << "\t" << sp.second << "\n";
            }
        }

        typedef Google< size_t, double >::Map probMapT;

        auto getOrElse = [] ( probMapT& pm, const size_t& key, double alt ) {
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
            auto alphas = computeAlphasDouble( countDict[*vit], k, total );

            // for each top-k score class
            for ( size_t i = 0; i < countDict[*vit].size(); ++i ) {
                double pScore;
                cl_I pCount;
                // The score and it's count
                tie(pScore, pCount) = countDict[*vit][i];

                double tprob = 0.0;
                // for all incoming edges contributing to this score
                for ( const auto& e : usedEdges[*vit][pScore] ) {
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
                    for ( const auto& tind : tail ) {
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
        for ( const auto& sc : countDict[rootInd] ) {
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


    vector<CostClass> computeKBest(const size_t& vid,
                                   const size_t& k,
                                   unordered_map<size_t, vector<CostClass> >& tkd,
                                   unique_ptr<ForwardHypergraph>& H) {

        // Dictionary that hold, for each incoming edge, the
        // number of score classes for each tail node
        unordered_map<size_t, vector<size_t>> sizeDict;


        // Priority queue of derivations for the given vertex
        using boost::heap::fibonacci_heap;
        QueueCmp<edvsT> ord;
        fibonacci_heap<edvsT, boost::heap::compare<QueueCmp<edvsT>>> vpq(ord);

        // Will hold the top-k cost classes for this solution
        vector<CostClass> cc;
        cc.reserve(k);

        auto round3 = [=]( double num ) -> double {
            double result = num * 1000;
            result = std::floor(result);
            result = result / 1000;
            return result;
        };

        for( auto& eid : H->incident(vid) ) {
            // The edge, backpointer array and cost of the derivation
            auto edge = H->edge(eid);
            vector<size_t> bp(edge.tail().size(), 0);
            auto cost = round3(getCost(tkd, bp, eid, H));

            // Push this potential derivation on the queue
            vpq.push( make_tuple( cost, eid, bp ) );

            // Fill in the size array for this edge
            vector<size_t> sizes;
            sizes.reserve(edge.tail().size());
            for( auto& tn : edge.tail() ) {
                sizes.push_back(tkd[tn].size());
            }
            sizeDict[eid] = sizes;
        }

        // Compute the cost of a derivation
        std::function< double( const size_t&, const vector<size_t>& ) > computeScore = [&] (const size_t& eid,
        const vector<size_t>& inds ) -> double {
            return round3( getCost(tkd, inds, eid, H) );
        };

        // Exact score classes
        double epsilon = 0.5;

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
                CostClass cclass(cost);
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
        if ( cc.size() == k+1 ) {
            cc.resize(k);
        }
        return cc;
    }

    void viterbiCountNew( unique_ptr<ForwardHypergraph>& H, TreePtrT& t, TreeInfo& ti, double penalty, const vector<size_t>& order,
                          slnDictT& slnDict, countDictT& countDict, const size_t& k,
                          const string& outputName, const vector<FlipKey>& outputKeys ) {

        // Compute the *weighted* probability of each edge being in
        // the top k distinct scoring solutions

        // Dictionary that holds the top-k cost classes for each
        // vertex, as well as other relevant information
        unordered_map< size_t, vector<CostClass> > tkd;

        auto vstr = [&]( const FlipKey& vert ) -> string {
            auto uname = t->getNodeName(vert.u());
            auto vname = t->getNodeName(vert.v());
            if (uname > vname) {
                auto tmp = vname;
                vname = uname;
                uname = tmp;
            }
            auto fstr = vert.f() ? "true" : "false";
            auto rstr = vert.r() ? "true" : "false";
            return uname+"\t"+vname+"\t"+fstr+"\t"+rstr;
        };

        double costsum = 0.0;

        // Each leaf has a single solution which is, by definition, of optimal cost
        for ( const auto& vit : order ) {
            if ( H->incident(vit).size() == 0 ) {
                tkd[vit] = { CostClass(slnDict[vit][0].cost) };
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
        size_t ctr=0;
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
            fs.set_empty_key( Derivation::flipT(-1,-1,"") );
            slnDict[*vit][0] = Derivation( tkd[*vit].front().cost(), 0, vector<size_t>(), fs);

        } // loop over verts

        typedef Google< size_t, double >::Map probMapT;

        probMapT probMap;
        probMapT outProbMap;
        probMap.set_empty_key( std::numeric_limits<size_t>::max() );
        outProbMap.set_empty_key( std::numeric_limits<size_t>::max() );

        // Map from a vertex to the maximum cost class that is
        // considered in deriving any solution that is actually used.
        vector<size_t> maxCostClass(order.size(),0);

        FlipKey rootKey( t->getRootId(), t->getRootId(), false, false);
        auto rootInd = H->index(rootKey);
        // We always consider all k classes for the root
        maxCostClass[rootInd] = k;
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
            auto parentProb = (probMap.find(*vit) == probMap.end() ? 0.0 : probMap[*vit];

            // The total # of derivations of this vertex (over all considered cost classes)
            cl_I total(0);
            vector< tuple<double, cl_I> > cd;
            //for ( size_t i = 0; i < tkd[*vit].size(); ++i) {
            for ( size_t i = 0; i < maxCostClass[*vit]; ++i ) {
                total += tkd[*vit][i].total();
                for( auto& e : tkd[*vit][i].usedEdges() ) {
                    auto frontier = tkd[*vit][i].getEdgeFrontier(e);

                    auto tail = H->edge(e).tail();
                    for( size_t j = 0; j < tail.size(); ++j) {
                        auto tn = tail[j];
                        maxCostClass[tn] = std::max( frontier[j]+1, maxCostClass[tn] );
                    }
                }
                cd.push_back( make_tuple(tkd[*vit][i].cost(), tkd[*vit][i].total()) );
            }

            // Compute the weights for each of the score classes
            // considered at this node
            auto alphas = computeAlphasDouble( cd, k, total );

            // for each top-k score class
            //for ( size_t i = 0; i < tkd[*vit].size(); ++i ) {
            for ( size_t i = 0; i < maxCostClass[*vit]; ++i ) {

                // The i-th cost class for this node
                auto cc = tkd[*vit][i];
                double tprob = 0.0;

                // for all incoming edges contributing to this cost class
                for ( const auto& e : cc.usedEdges() ) {

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
                    for ( const auto& tind : tail ) {
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
        for ( const auto& cc : tkd[rootInd] ) {
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
                    cout << "inProbMap has " << probMap[vid] << ", outProbMap has" << outProbMap[vid] << "\n";
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

    vector<cl_RA> computeAlphas( const vector<tuple<double, cl_I>>& slnVec ) {
        vector<cl_RA> scores;
        cl_RA invSum(0);
        cl_RA one(1);
        for ( const auto& e : slnVec) {
            scores.push_back( one / static_cast<size_t>(get<0>(e)+1.0) );
            invSum += scores.back();
        }
        //cl_RA invSum = 1 / sum;
        // cerr << "invSum = " << invSum << "\n";
        vector<cl_RA> alphas;
        alphas.reserve(slnVec.size());
        for ( const auto& s : scores) {
            alphas.push_back( s/invSum );
        }
        return alphas;
    }


    double estimateThermodynamicBeta( const vector<tuple<double, cl_I>>& slnVecIn, const double& emin ) {
        if (slnVecIn.size() == 1) {
            return 0.0;
        }
        vector< tuple<double, mpreal> > slnVec;
        slnVec.reserve(slnVecIn.size());
        std::stringstream ss;
        for ( const auto& e : slnVecIn ) {
            ss << get<1>(e);
            mpreal n = ss.str();
            ss.clear();
            slnVec.push_back( make_tuple( get<0>(e), n ) );
        }
        double ctr = 0.0;
        double beta = 0.0;
        size_t n = slnVec.size();
        double scale = 1.0 / ( n-1 );
        for ( auto i = slnVec.begin(); (i+1) != slnVec.end(); ++i ) {
            auto j = i+1;
            auto lnI = mpfr::log( get<1>(*i) );
            auto eI = get<0>(*i) - emin;
            auto lnJ = mpfr::log( get<1>(*j) );
            auto eJ = get<0>(*j) - emin;
            double denom = (eJ-eI);

            if ( denom >= 1.0 ) {
                beta += scale * ( (lnJ - lnI).toDouble() ) / denom;
                ctr += 1.0;
            }
        }

        //beta /= ctr;
        //std::cerr << " ===== beta = " << beta << " ===== \n";
        return beta;
    }

    vector<double> getAlphas( const vector< tuple<double, cl_I> >& slnVec, const cl_I& numPaths ) {
        typedef tuple<double,cl_I> elemT;

        auto minmax = std::minmax_element( slnVec.begin(), slnVec.end(), [] (const elemT& e0, const elemT& e1) {
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
        for ( const auto& e : slnVec) {
            invScores.push_back( (alpha*double_approx(get<1>(e) / numPaths)) + ((1.0 - alpha) * std::exp( bestScore-get<0>(e) * scale ))  );
            totalWeight += invScores.back();
        }

        vector<double> alphas;
        alphas.reserve(slnVec.size());
        for ( const auto& e : invScores ) {
            alphas.push_back( e / totalWeight );
        }
        return alphas;
    }



    void insideOutside( unique_ptr<ForwardHypergraph>& H, TreePtrT& t, TreeInfo& ti, double penalty, const vector<size_t>& order, slnDictT& slnDict, countDictT& countDict, const string& outputName ) {

        // We'll use vectors for these since the nodes and vectors have a well
        // defined ordering.
        vector< tuple<double,cl_I> > insideScores( H->order() );
        vector< tuple<double,cl_I> > edgeWeightMap( H->size() );

        // Each leaf has a single solution which is, by definition, of optimal cost
        for ( const auto& vit : order ) {
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
                for ( const auto& e : H->incident(*vit) ) {

                    auto tail = H->getTail(e);
                    vector< tuple<double, cl_I> > tailScores;
                    tailScores.reserve( tail.size() );

                    cl_I numPaths(1);
                    cl_I summedNumPaths(0);
                    // Loop over all tail vertices
                    for ( const auto& tind : tail ) {
                        tailScores.push_back( insideScores[tind] );
                        numPaths *= get<1>(insideScores[tind]);
                        summedNumPaths += get<1>(insideScores[tind]);
                        cerr << "tailNode # paths = " << get<1>(insideScores[tind]) << "\n";
                    }
                    cerr << "summedNumPaths = " << summedNumPaths << "\n";
                    totalPaths += numPaths;
                    double avgTailCost = 0.0;
                    size_t tind = 0;
                    vector<double> talphas(getAlphas(tailScores,summedNumPaths));
                    for ( const auto& s : tailScores ) {
                        //cerr << "(score, frac) = " << get<0>(s) << ", (" << talphas[tind] << ")\n";
                        avgTailCost += get<0>(s);// * (1.0 / tailScores.size());//talphas[tind];
                        tind += 1;
                    }
                    //if (std::isnan(avgTailCost) || ctr > 200 ) { std::abort(); }
                    scores.push_back( make_tuple( H->edge(e).weight() + avgTailCost , numPaths ) );
                }

                double avgHeadScore = 0.0;
                for ( const auto& s : scores ) {
                    avgHeadScore += get<0>(s) * 1.0 / scores.size();//double_approx( get<1>(s) / totalPaths );
                }

                vector<double> contribWeights(getAlphas(scores,totalPaths));
                insideScores[*vit] = make_tuple(avgHeadScore, totalPaths);

                // once we know the score of the parent node, we can
                // compute the relative contribution from each incoming edge
                size_t eind = 0;
                for ( const auto& e : H->incident(*vit) ) {
                    edgeWeightMap[ e ] = make_tuple(contribWeights[eind], get<1>(scores[eind]));
                    ++eind;
                }
            }
        }

        typedef Google< size_t, double >::Map probMapT;

        auto getOrElse = [] ( probMapT& pm, const size_t& key, double alt ) {
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

            for ( const auto& e : H->incident(*vit) ) {
                auto condProb = get<0>(edgeWeightMap[e]);
                // for all tail vertices of this edge
                auto tail = H->getTail(e);
                for ( const auto& tind : tail ) {
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

    FlipKey keyForAction( const FlipKey& fk , const string& ft ) {
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





    void viterbi( unique_ptr<ForwardHypergraph>& H, TreePtrT& t, TreeInfo& ti, double penalty, const vector<size_t>& order, slnDictT& slnDict ) {
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

for ( auto& e : H->incident(*vit) ) {
                auto edge = H->edge(e);
                auto hind = edge.head();
                auto w = edge.weight();
                auto tvert = H->vertex(edge.tail()[0]);

                Google<Derivation::flipT>::Set flips;
                flips.set_empty_key( make_tuple(-1,-1,""));
                //set<Derivation::flipT> flips;
                string children;
for ( auto tn : edge.tail() ) {
                    auto cvert = H->vertex(tn);
                    w += slnDict[tn][0].cost;
                    /*std::set_union( slnDict[tn][0].flips.begin(), slnDict[tn][0].flips.end(),
                                    slnDict[tn][0].flips.begin(), slnDict[tn][0].flips.end(),
                                    std::inserter( flips, flips.begin() ) );
                    */
for ( auto& f : slnDict[tn][0].flips ) {
                        flips.insert(f);
                    }

                    //children += "[ " + t->getNodeName(cvert.u()) + ", " + t->getNodeName(cvert.v()) + " : (" + lexical_cast<string>(cvert.f()) + lexical_cast<string>(cvert.r()) +")] ";
                }

                auto ft = flipType(vert,tvert);
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
for ( auto& e : v ) {
            cerr << e << " ";
        }
        cerr << "]";
    }

    void getCandidates( const unique_ptr<ForwardHypergraph>& H, TreeInfo& ti, double penalty, slnDictT& d, size_t v, size_t k, vector<Derivation>* kbest ) {
        // There is only one candidate for any leaf node
        auto incoming = H->incident(v);
        auto numEdges = incoming.size();
        if ( numEdges == 0 ) {
            return;
        }
        // For internal nodes
        // List of potential candidates
        vector<Derivation> temp;
        // The head vertex and it's constituent tree nodes
        auto hvert = H->vertex(v);
        auto a = hvert.u();
        auto b = hvert.v();

        // Initalize the heap with the 1-best solution for each incoming edge
for ( auto e : incoming ) {
            // Compute the score of this derivation
            auto edge = H->edge(e);
            auto w = edge.weight();
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

for ( auto t : tail ) {
                s += d[t][0].cost;
                /*std::set_union( d[t][0].flips.begin(), d[t][0].flips.end(),
                                d[t][0].flips.begin(), d[t][0].flips.end(),
                                std::inserter( ed, ed.begin() ) );
                */
for ( auto& f : d[t][0].flips ) {
                    ed.insert(f);
                }
            }

            if ( ft != "n" ) {
                ed.insert( make_tuple(a,b,ft) );
            }

            //make_heap( temp.begin(), temp.end() );
            temp.push_back( Derivation(s, e, j, ed) );
            //push_heap( temp.begin(), temp.end() );

        }
        // We only need the top k solutions here
        // kbest = heapq.nsmallest( max(k,numEdges), temp );
        // Heapify the list and return it
        k = std::min( temp.size(), std::max( k, numEdges ) );
        partial_sort( temp.begin(), temp.begin()+k, temp.end(),
        [] ( const Derivation& d1, const Derivation& d2 ) {
            return d1.cost < d2.cost;
        } );
        kbest->reserve(k);
        copy( temp.begin(), temp.begin()+k, std::back_inserter(*kbest) );
        //for ( auto c = temp.begin(); c != temp.begin()+k; ++c ) {
        //    kbest.push_back(*c);
        //}
        //cerr << "|temp| = " << temp.size() << ", |kbest| = " << kbest.size() << "\n";
    }

    bool sameEffect( const Derivation& d0, const Derivation& d1 ) {
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
        [&] ( const tuple<int,int,string>& d ) {
            return d1.flips.find(d) == d1.flips.end();
        } );

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
        auto a = hvert.u();
        auto b = hvert.v();

        auto ft = flipType(hvert, tvert);

        // For each tail node in this edge
        size_t i = 0;
for ( auto vind : tail ) {
            //cerr << "before : deriv ";
            //printVector(j);
            //cerr << "\n";

            // Ask for the next best derivation at this tail node
            vector<size_t> jp(j);
            jp[i] += 1;

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
                auto s = w;
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
                    auto vp = tail[ind];
                    auto k = jp[ind];
                    s += d[ vp ][ k ].cost;
for ( auto& f : d[vp][k].flips ) {
                        ed.insert(f);
                    }
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
                    if ( sameEffect(deriv, od) ) {
                        foundDer = true;
                        break;
                    }
                }
                if ( !foundDer ) {
                    locCand->push_back( deriv );
                    push_heap( locCand->begin(), locCand->end(),
                    [] ( const Derivation& d0, const Derivation& d1) {
                        return d0.cost > d1.cost;
                    } );
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
                [] ( const TaggedDerivT& a, const TaggedDerivT& b ) {
                    return a.first < b.first;
                } )->first;
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
                [] ( const Derivation& d0, const Derivation& d1) {
                    return d0.cost > d1.cost;
                } );
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
                auto existingElem = find_if( d[v].begin(), d[v].end(), [&] ( const TaggedDerivT& od ) {
                    return sameEffect(deriv,od.second);
                } );
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

using GraphUtils::undirectedGraphT;
using GraphUtils::directedGraphT;

template void MultiOpt::leafCostDict< undirectedGraphT  >( std::unique_ptr<ForwardHypergraph>& , Utils::Trees::TreePtrT& , undirectedGraphT& , bool , double , double , MultiOpt::slnDictT& );

template void MultiOpt::leafCostDict< directedGraphT >( std::unique_ptr<ForwardHypergraph>& , Utils::Trees::TreePtrT& , directedGraphT& , bool , double , double , MultiOpt::slnDictT& );
