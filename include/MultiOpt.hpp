#ifndef MULTI_OPT_HPP
#define MULTI_OPT_HPP

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
#include <limits>

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

    using google::dense_hash_set;
    using google::dense_hash_map;
    using google::sparse_hash_map;
    using std::unordered_map;
    using std::unordered_set;
    using std::tuple;
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

    typedef unordered_map< size_t, unordered_map<size_t, Derivation>> slnDictT;
    //typedef unordered_map<size_t, vector<Derivation>> slnDictT;
    //typedef Google<size_t, vector<Derivation>>::Map slnDictT;
    //typedef Google<size_t, Google<size_t,Derivation>::Map>::Map slnDictT;

    auto none = make_tuple(false,false); auto fw = make_tuple(true,false); auto rev = make_tuple(false,true); auto both = make_tuple(true,true);

    static FlipMapT flipDict = {
        {none, { {both, "b+"}, {fw, "f+"}, {rev, "r+"}, {none, "n"} } },
        {fw,   { {both, "r+"}, {rev, "f-r+"}, {none, "f-"}, {fw, "n"} } },
        {rev,  { {both, "f+"}, {rev, "n"}, {none, "r-"}, {fw, "f+r-"} } },
        {both, { {both, "n"}, {fw, "r-"}, {rev, "f-"}, {none, "b-"} } }
    };

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

    unique_ptr<ForwardHypergraph>  buildSolutionSpaceGraph( const Tree* t, const TreeInfo& ti, bool directed ) {
        Google<int>::Set leafSet;
        leafSet.set_empty_key(-1);
        for ( auto l : t->getLeavesId() ){ leafSet.insert(l); }
        auto isLeaf = [&]( const int& nid ) { return leafSet.find(nid) != leafSet.end(); };//t->isLeaf(nid); };
        auto isInternal = [&]( const int& nid ) { return ! isLeaf(nid); };
        auto isLost = [&]( int nid ){ return (t->getNodeName(nid)).find("LOST") != std::string::npos; };

        unique_ptr<ForwardHypergraph> slnSpaceGraph( new ForwardHypergraph() );

        auto addIncomingHyperedge = [=,&slnSpaceGraph] (const FlipKey& k, const int& rnode, const int& onode ) {
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
                    if ( dualFlipEdge.size() > 0 ) { slnSpaceGraph->addEdge( dualFlipEdge, k, 1.0 ); }
                }
            } else {
                if ( isInternal(rnode) ) {
                    //cout << t->getNodeName(rnode) << " is INTERNAL\n";
                    // Get the children nodes
                    int LRN = t->getSonsId(rnode)[0]; int RRN = t->getSonsId(rnode)[1];
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
                        slnSpaceGraph->addEdge( noFlip, k, 0.0 );
                    }
                    if ( dualFlip.size() > 0 ) {
                        auto w = directed ? 2.0 : 1.0;//2.0 ? directed : 1.0;
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
                            slnSpaceGraph->addEdge( fwFlip, k , 1.0 );
                        }
                        if ( revFlip.size() > 0 ) {
                            slnSpaceGraph->addEdge( revFlip, k, 1.0 );
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
        //cout << t->getNumberOfNodes() << "#";
        //cout << nodes.size() << "\n";
        auto tbegin = nodes.cbegin(); auto tend = nodes.cend();
        vector<int>::const_iterator uit = tbegin; vector<int>::const_iterator vit = tbegin;

        //for( auto l : nodes ) { cout << l << " ";}
        //cout << "\n\n";

        for ( uit = tbegin; uit != tend; uit++ ) {
            int u = *uit;
            for ( vit = uit; vit != tend; vit++ ) {
                int v = *vit;
                //cout << "u,v = " << u << ", " << v << "\n";
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

        //cout << "Added vertices\n";

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
        return slnSpaceGraph;
    }



    template< typename GT >
    void leafCostDict( unique_ptr<ForwardHypergraph>& H, Tree* T, GT& G, bool directed, double cc, double dc, slnDictT& slnDict ) {
        //cout << "In leafCostDict\n";
        /*
          Given the duplication tree T, the root vertex rv, the extant graph G and
          the constraints, fill in the values for the leaf nodes of the hypergraph
        */
        auto undirected = !directed;
        // Cost of going from the state encoded by a node to the state of a true graph edge
        typedef tuple<bool,bool> flipTupleT;
        typedef tuple<double, string> costRepT;
        typedef unordered_map< flipTupleT, unordered_map<flipTupleT, costRepT> > costMapT;
        typedef unordered_map< flipTupleT, unordered_map<bool, costRepT> > selfCostMapT;

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

        selfCostMapT selfLoopCostDict;
        selfLoopCostDict[ none ][ true ] = make_tuple( cc, "b+" );
        selfLoopCostDict[ none ][ false ] = make_tuple( 0.0, "n" );

        selfLoopCostDict[ rev ][ true ] = make_tuple( 0.0, "n" );
        selfLoopCostDict[ rev ][ false ] = make_tuple( dc, "b-" );

        selfLoopCostDict[ fw ][ true ] = make_tuple( 0.0, "n" );
        selfLoopCostDict[ fw ][ false ] = make_tuple( dc, "b-" );

        selfLoopCostDict[ both ][ true ] = make_tuple( 0.0, "n" );
        selfLoopCostDict[ both ][ false ] = make_tuple( dc, "b-" );

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
        /*
          for ( auto lhn : leafHypernodes ) {
          auto ve = H->vertex(lhn);
          auto u = T->getNodeName(ve.u()); auto v = T->getNodeName(ve.v());
          auto f = ve.f(); auto r = ve.r();
          cout << u << ", " << v << " :: [(" << f << ", " << r << ")]\n";
          } cout << "\n";
        */

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
        //auto leaves = T->getLeavesId();
        for ( auto l : T->getLeavesId() ) {
            leafIds.insert(l);
            //cout << "leaf " << l << "\n";
        }

        auto vp = boost::vertices(G);
        for ( auto it = vp.first; it != vp.second; ++it ) {
            auto v = *it;
            auto idx = G[v].idx;
            // found this node's id in the set of extant vertices'
            if ( leafIds.find(idx) != leafIds.end() ) { extantNodes.insert(v); }
        }
        /*
          cout << "extant nodes : \n";
          for ( auto n : extantNodes ) { cout << T->getNodeName( G[n].idx ) << " ";} cout << "\n";
        */

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
                // set<tuple<int,int,string>> es;
                slnDict[n] = { {0,Derivation(lostCost, n, ev, es)} };
            } else {
                // Otherwise, u and v both exist in the extant
                // network, so get the appropriate info
                auto u = idToVertMap[nd.u()];
                auto v = idToVertMap[nd.v()];
                auto f = nd.f(); auto r = nd.r();
                //cout << " HN [" << T->getNodeName(u) << ", " << T->getNodeName(v) << ", (" << f << r << ")]\n";
                if (u != v) {
                    auto d_f = hasEdge(u,v);
                    auto d_r = hasEdge(v,u);

                    auto costFlip = costDict[ make_tuple(f,r) ][ make_tuple(d_f,d_r) ];
                    auto cost = get<0>(costFlip); auto flip = get<1>(costFlip);
                    //cout << "cost = " << cost << ", flip = " <<
                    //flip;
                    Google<Derivation::flipT>::Set effectiveEdges;
                    effectiveEdges.set_empty_key( make_tuple(-1,-1,""));
                    //set< tuple<int,int,string> > effectiveEdges;
                    if ( flip != "n" ) { effectiveEdges.insert( make_tuple(nd.u(),nd.v(),flip) ); }
                    vector<size_t> ev;
                    slnDict[n] = { {0,Derivation(cost, n, ev, effectiveEdges)} };

                } else {
                    auto hasSelfLoop = hasEdge(u,v);
                    auto costFlip = selfLoopCostDict[ make_tuple(f,r) ][ hasSelfLoop ];
                    auto cost = get<0>(costFlip); auto flip = get<1>(costFlip);
                    Google<Derivation::flipT>::Set effectiveEdges;
                    effectiveEdges.set_empty_key( make_tuple(-1,-1,""));
                    //set< tuple<int,int,string> > effectiveEdges;
                    if ( flip != "n" ) { effectiveEdges.insert( make_tuple(nd.u(),nd.v(),flip) ); }
                    vector<size_t> ev;
                    slnDict[n] = { {0,Derivation(cost, n, ev, effectiveEdges)} };

                } // ( u != v )
            } // ( lostU || lostV )
        } // loop over leaf hypernodes
    }


    void viterbi( unique_ptr<ForwardHypergraph>& H, Tree* t, const vector<size_t>& order, slnDictT& slnDict ) {
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

                if ( w > 0 ) { flips.insert( make_tuple(vert.u(), vert.v(), flipType(vert,tvert)) ); }
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

    void getCandidates( const unique_ptr<ForwardHypergraph>& H, slnDictT& d, size_t v, size_t k, vector<Derivation>* kbest ){
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

            for ( auto t : tail ){
                s += d[t][0].cost;
                /*std::set_union( d[t][0].flips.begin(), d[t][0].flips.end(),
                                d[t][0].flips.begin(), d[t][0].flips.end(),
                                std::inserter( ed, ed.begin() ) );
                */
                for ( auto& f : d[t][0].flips ) { ed.insert(f); }
            }

            if ( w > 0.0 ) { ed.insert( make_tuple(a,b,ft) ); }

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

    void lazyKthBest( unique_ptr<ForwardHypergraph>& H, slnDictT& d, size_t v, size_t k, size_t kp );

    void lazyNext( unique_ptr<ForwardHypergraph>& H, slnDictT& d, vector<Derivation>* locCand, size_t eind,
                   const vector<size_t>& j, size_t kp) {

        auto edge = H->edge(eind);
        auto tail = edge.tail();
        auto hvert = H->vertex(edge.head());
        auto tvert = H->vertex(tail[0]);
        auto a = hvert.u(); auto b = hvert.v();

        auto ft = flipType(hvert, tvert);

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
            lazyKthBest( H, d, vind, jp[i], kp );



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
                    auto vp = tail[ind]; auto k = jp[ind];
                    s += d[ vp ][ k ].cost;
                    for ( auto& f : d[vp][k].flips ) { ed.insert(f); }
                }


                if ( w > 0.0 ) {
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


    void lazyKthBest( unique_ptr<ForwardHypergraph>& H, slnDictT& d, size_t v, size_t k, size_t kp ) {
        // If this is our first vist to vertex v, then
        // populate its candidate list

        //cerr << "Searching for candidates for " << v << "\n";
        if ( cand.find(v) == cand.end() ) {
            //cerr << "There were none; adding to the cand map \n";
            //if ( kp == 0 ) { recursedStore.set_empty_key( make_tuple(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()) ); }
            cand[v] = new vector< Derivation >();
            //cerr << "done \n";
            getCandidates(H, d, v, kp, cand[v]);
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
                        lazyNext(H, d, cand[v], deriv.target, deriv.bp, kp);
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
                        lazyNext(H, d, cand[v], deriv.target, deriv.bp, kp);
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
