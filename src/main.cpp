/** Standard Includes */
#include <map>
#include <set>
#include <stdexcept>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <fstream>

/** Boost Includes */
#include <boost/program_options.hpp>
#include <boost/progress.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

/** Local Includes */
#include <arg_parser.hpp>

/** Bio++ **/
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Exceptions.h>

/** Boost **/
#include <algorithm>                 // for std::for_each
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/variant.hpp>
#include "GraphUtils.hpp"
#include "FlipKey.hpp"
#include "MultiOpt.hpp"
#include "utils.hpp"

/** Namespace uses */
using std::string;
using std::ifstream;
using std::map;
using std::ofstream;
using std::ios;
using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::exception;
using std::abort;
using std::unordered_map;
using std::unordered_set;
using std::set;
using bpp::Newick;
using bpp::Tree;
using bpp::Exception;
using boost::adjacency_list;
using boost::graph_traits;
using boost::vecS;
using boost::listS;
using boost::variant;
using Utils::TreeInfo;

const std::string getName(Tree* t, int nid) {
    if (t->hasNodeName(nid)) {
        return t->getNodeName(nid);
    } else {
        return reinterpret_cast<bpp::BppString*>( t->getNodeProperty(nid,"name"))->toSTL();
    }
}

string getExtantNetwork(const string& s) {
    return ( s.find("LOST") != string::npos ) ? "LOST" : s.substr( s.find_last_of('_') );
};


void prepareTree( Tree* t, TreeInfo& ti, int nid ) {
    if (t->isLeaf(nid)) {
        ti.leaves[nid] = {nid};
        ti.subnodes[nid] = {nid};
        ti.enets[nid] = {getExtantNetwork(getName(t,nid))};

    } else {
        ti.leaves[nid] = { nid };
        ti.subnodes[nid] = { nid };
        ti.enets[nid] = {};

        for ( auto cid : t->getSonsId(nid) ) {
            prepareTree( t, ti, cid );
        }

        for ( auto cid : t->getSonsId(nid) ) {
            for ( auto l : ti.leaves[cid] ) { ti.leaves[nid].insert(l); }
            for ( auto l : ti.subnodes[cid] ) { ti.subnodes[nid].insert(l); }
            for ( auto l : ti.enets[cid] ) { ti.enets[nid].insert(l); }
        }
    }
    /*
    cout << "subnodes of " << nid << "\n";
    for( auto n : ti.subnodes[nid] ) {
        cout << n << " ";
    }cout << "\n";
    */
}


int main( int argc, char **argv ){
    po::options_description desc("Allowed Options");

    desc.add_options()
        ("help,h", "produce this message")
        ("dupHist,d", po::value< string >(), "mesh input file")
        ("target,t", po::value< string >(), "extant graph file")
        ("undir,u", po::value< bool >()->zero_tokens() , "graph is undirected")
        ;

    try {
        ArgParser ap( desc );
        ap.parse_args( argc, argv );
        string treeName = ap["dupHist"].as<string>();
        string graphName = ap["target"].as<string>();

        bool undirected = ap.isPresent("undir");
        bool directed = !undirected;
        cout << "UNDIRECTED = " << undirected << "\n";

        auto newickReader = new Newick(true,true); //No comment allowed!
        newickReader->enableExtendedBootstrapProperty("name");
        try {
            auto tree = newickReader->read(treeName); // Tree in file
            // MyTestTree.dnd

            cout << "Tree has " << tree->getNumberOfNodes() << " nodes." << endl;
            //auto rootId = tree->getRootId();

            if ( tree->isRooted() ){ cout << "Tree is rooted\n"; }
            for ( auto nid : tree->getNodesId() ) {
                for ( auto prop : tree->getBranchPropertyNames(nid) ) {
                    tree->setNodeName(nid, reinterpret_cast<bpp::BppString*>(tree->getBranchProperty(nid,prop))->toSTL() );
                }
            }
            /*
            for (auto nid : tree->getNodesId()) {
                if (tree->hasNodeName(nid)) {
                    cout << "Node : " << tree->getNodeName(nid) << "\n";
                } else {
                    for (auto prop : tree->getNodePropertyNames(nid)) {
                        auto p = reinterpret_cast<bpp::BppString*>(tree->getNodeProperty(nid,prop));
                        cout << "Node : " << p->toSTL() << "\n";
                    }
                }
            }
            */
            TreeInfo tinfo("TreeInfo");
            prepareTree( tree, tinfo, tree->getRootId() );

            for ( auto nid : tree->getNodesId() ){
                auto nodeName = getName(tree,nid);
                cout << "Node [" << nodeName << "]\n";
                for ( auto enet : tinfo.enets[nid] ) {
                    cout << enet << ", ";
                }
            }


            unordered_map<string, int> leafIDMap;
            for( auto nid : tree->getLeavesId() ) { leafIDMap[ getName(tree,nid) ] = nid; }

            typedef adjacency_list<boost::hash_setS, vecS, boost::directedS, GraphUtils::Node > directedGraphT;
            typedef adjacency_list<boost::hash_setS, vecS, boost::undirectedS, GraphUtils::Node > undirectedGraphT;
            variant< directedGraphT, undirectedGraphT > G;
            if ( undirected ) {
                G = undirectedGraphT();
                GraphUtils::readFromAdjacencyList(graphName, leafIDMap, get<undirectedGraphT>(G) );
                auto vp = vertices(get<undirectedGraphT>(G));
                size_t i = 0;
                for( auto it = vp.first; it != vp.second; ++it ) { cout << "Vertex " << *it << ", name =  " << get<undirectedGraphT>(G)[*it].name << " # " << i << "\n"; ++i; }
            } else {
                G = directedGraphT();
                GraphUtils::readFromAdjacencyList(graphName, leafIDMap, get<directedGraphT>(G) );
            }

            //typedef graph_traits < graphT >::vertex_descriptor vertexDescT;
            //typedef graph_traits < graphT >::edge_descriptor edgeDescT;

            auto H = MultiOpt::buildSolutionSpaceGraph( tree, tinfo, directed );
            vector<size_t> order; order.reserve( H->order() );
            MultiOpt::topologicalOrder( H, order );

            MultiOpt::slnDictT slnDict;
            if ( undirected ) {
                MultiOpt::leafCostDict( H, tree, get<undirectedGraphT>(G), directed, 1.0, 1.0, slnDict);
            } else {
                MultiOpt::leafCostDict( H, tree, get<directedGraphT>(G), directed, 1.0, 1.0, slnDict);
            }

            MultiOpt::viterbi( H, tree, order, slnDict );

            auto rootKey = FlipKey( tree->getRootId(), tree->getRootId(), false, false );
            auto rootInd = H->index(rootKey);

            auto totalDerivs = 1000;
            auto derivNum = 0;
            double prevCost, newCost;
            prevCost = newCost = slnDict[ rootInd ][ derivNum ].cost;
            MultiOpt::initKBest();

            while ( derivNum < totalDerivs ) { // || prevCost == newCost ) {
                MultiOpt::lazyKthBest( H, slnDict, rootInd, derivNum, derivNum );
                cerr << "returned from call to lazy" << derivNum << "best!!\n";
                auto d = slnDict[ rootInd ][ derivNum ];
                cout << "DERIVATION # " << derivNum << ": for [(" << tree->getNodeName(rootKey.u()) << ", " << tree->getNodeName(rootKey.v()) <<  ") : (" << rootKey.f() << rootKey.r() << ")] (" << d.cost << ")\n";

                std::fstream output( "output.txt", std::fstream::in | std::fstream::out | std::fstream::app );
                for ( auto f : d.flips ) {
                    cout << " [ " << tree->getNodeName(get<0>(f)) << ", " << tree->getNodeName(get<1>(f)) << " : " << get<2>(f) <<  "]\n";
                    output << tree->getNodeName(get<0>(f)) << "\t" << tree->getNodeName(get<1>(f)) << "\t" << get<2>(f).substr(0,1) << "\n";
                }
                output.close();

                prevCost = newCost;
                newCost = d.cost;

                derivNum++;
            }

            delete tree;
        } catch (Exception e) {
            cout << "Error when reading tree." << endl;
        }
        delete newickReader;

    } catch(exception& e) {

        cerr << "Caught Exception: [" << e.what() << "]\n";
        abort();

    }

    return 0;
}
