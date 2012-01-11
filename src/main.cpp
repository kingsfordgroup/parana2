/** Standard Includes */
#include <map>
#include <set>
#include <stdexcept>
#include <cmath>
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
#include <boost/lexical_cast.hpp>

/** Local Includes */
#include <arg_parser.hpp>

/** Bio++ **/
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Io/NexusIOTree.h>
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

/** CLN **/
#include <cln/io.h>
#include <cln/rational_io.h>
#include <cln/integer_io.h>

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
using bpp::NexusIOTree;
using bpp::Tree;
using bpp::Exception;
using boost::adjacency_list;
using boost::graph_traits;
using boost::vecS;
using boost::listS;
using boost::variant;
using boost::lexical_cast;
using Utils::TreeInfo;

#include <cstdlib>

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
    // if the current node is not the root
    if ( nid != t->getRootId() ) {
        auto fid = t->getFatherId(nid);
        auto parentDeathT = get<1>(ti.extantInterval[fid]);
        assert ( ti.extantInterval.find(fid) != ti.extantInterval.end() );
        double fdist = 0.0;
        if ( t->isLeaf(nid) ) {
            fdist = std::numeric_limits<double>::infinity();
        } else if ( t->hasDistanceToFather(nid) ) {
            fdist = t->getDistanceToFather(nid);
        }
        ti.extantInterval[ nid ] = make_tuple( parentDeathT, parentDeathT + fdist );
    }


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
        ("ratio,r", po::value< double >()->default_value(1.0), "ratio of creation to deletion cost")
        ("undir,u", po::value< bool >()->zero_tokens() , "graph is undirected")
        ("output,o", po::value< string >()->default_value("edgeProbs.txt"), "output file containing edge probabilities")
        ("numOpt,k", po::value< size_t>()->default_value(10), "number of near-optimal score classes to use")
        ("timePenalty,p", po::value< double >()->default_value(0.0), "amount to penalize flips between nodes whose time intervals don't overlap'")
        ;

    try {
        ArgParser ap( desc );
        ap.parse_args( argc, argv );
        string treeName = ap["dupHist"].as<string>();
        string graphName = ap["target"].as<string>();
        string outputName = ap["output"].as<string>();

        double penalty = ap["timePenalty"].as<double>();
        double creationCost = ap["ratio"].as<double>();
        double deletionCost = 1.0;

        bool undirected = ap.isPresent("undir");
        bool directed = !undirected;
        size_t k = ap["numOpt"].as<size_t>();

        cout << "UNDIRECTED = " << undirected << "\n";

        auto newickReader = new Newick(true,true); //No comment allowed!
        newickReader->enableExtendedBootstrapProperty("name");
        try {
            auto nexusReader = new NexusIOTree();
            bpp::TreeTemplate<bpp::Node>* tree;
            try {
                tree = newickReader->read(treeName); // Tree in file
            } catch ( const exception &e ) {
                cerr << "ERROR READING TREE " << e.what();
            }
            cout << "Tree has " << tree->getNumberOfNodes() << " nodes." << endl;
            //auto rootId = tree->getRootId();

            if ( tree->isRooted() ){ cout << "Tree is rooted\n"; }
            for ( auto nid : tree->getNodesId() ) {
                for ( auto prop : tree->getBranchPropertyNames(nid) ) {
                    tree->setNodeName(nid, reinterpret_cast<bpp::BppString*>(tree->getBranchProperty(nid,prop))->toSTL() );
                }
            }

            // cout how many branches had length 0
            /*
            int zcount = 0;
            for ( auto n : tree->getNodesId() ) {
                if ( n != tree->getRootId() ) {
                    auto dist = tree->getDistanceToFather(n);
                    if ( dist == 0.0 ) { zcount++; }
                    cout << tree->getDistanceToFather(n) << " ";
                }
            } cout << "\n";
            cout << zcount << " branches had length 0\n";
            */
            TreeInfo tinfo("TreeInfo");
            auto rId = tree->getRootId();
            tinfo.extantInterval[ rId ] = make_tuple( -std::numeric_limits<double>::infinity(), 0.0 );
            prepareTree( tree, tinfo, rId );

            // print out the existence intervals for all of the nodes

            for ( auto nid : tree->getNodesId() ){
                auto nodeName = getName(tree,nid);
                cout << "Node [" << nodeName << "] : existance interval = [" <<
                    get<0>(tinfo.extantInterval[nid]) << ", " << get<1>(tinfo.extantInterval[nid]) << "] \n";
            }



            unordered_map<string, int> leafIDMap;
            for( auto nid : tree->getLeavesId() ) { leafIDMap[ getName(tree,nid) ] = nid; }

            auto isAdjList = [&](){
                auto spos = graphName.find_last_of(".");
                auto suffix(graphName.substr( spos ));
                return suffix == ".adj";
            }();

            typedef adjacency_list<boost::hash_setS, vecS, boost::directedS, GraphUtils::Node, GraphUtils::Edge > directedGraphT;
            typedef adjacency_list<boost::hash_setS, vecS, boost::undirectedS, GraphUtils::Node, GraphUtils::Edge > undirectedGraphT;
            variant< directedGraphT, undirectedGraphT > G;
            if ( undirected ) {
                G = undirectedGraphT();
                if ( isAdjList ) {
                    GraphUtils::readFromAdjacencyList( graphName, leafIDMap, get<undirectedGraphT>(G) );
                } else {
                    GraphUtils::readFromMultilineAdjacencyList(graphName, leafIDMap, get<undirectedGraphT>(G) );
                }

                /*
                auto vp = vertices(get<undirectedGraphT>(G));
                size_t i = 0;
                for( auto it = vp.first; it != vp.second; ++it ) {cout << "Vertex " << *it << ", name =  " << get<undirectedGraphT>(G)[*it].name << " # " << i << "\n"; ++i; }
                */
            } else {
                G = directedGraphT();
                GraphUtils::readFromMultilineAdjacencyList(graphName, leafIDMap, get<directedGraphT>(G) );
            }

            //typedef graph_traits < graphT >::vertex_descriptor vertexDescT;
            //typedef graph_traits < graphT >::edge_descriptor
            //edgeDescT;

            auto H = MultiOpt::buildSolutionSpaceGraph( tree, tinfo, creationCost, deletionCost, penalty, directed );

            vector<size_t> order; order.reserve( H->order() );
            MultiOpt::topologicalOrder( H, order );

            MultiOpt::slnDictT slnDict;
            //slnDict.set_empty_key( std::numeric_limits<size_t>::max() );
            if ( undirected ) {
                MultiOpt::leafCostDict( H, tree, get<undirectedGraphT>(G), directed, creationCost, deletionCost, slnDict);
            } else {
                MultiOpt::leafCostDict( H, tree, get<directedGraphT>(G), directed, creationCost, deletionCost, slnDict);
            }


            auto rootKey = FlipKey( tree->getRootId(), tree->getRootId(), false, false );
            auto rootInd = H->index(rootKey);

            // Count the # of opt slns.
            MultiOpt::countDictT countDict;
            countDict.set_empty_key(-1);
            MultiOpt::viterbiCount(H, tree, tinfo, penalty, order, slnDict, countDict, k, outputName);

            //cout << "The optimal cost solutions have a cost of " << slnDict[rootInd][0].cost << "\n";
            //cout << "There are " << get<1>(countDict[rootInd][0]) << " optimal solutions ";
            //cout << "Near optimal -- (cost, count) : \n";
            //for ( const auto& der : countDict[rootInd] ) {
            //    cout << " (" << get<0>(der) << ",  " << get<1>(der) << ")\n";
            //}
            //cout << "\n";
            return 0;


            MultiOpt::viterbi( H, tree, tinfo, penalty, order, slnDict );

            auto totalDerivs = 1;
            auto derivNum = 0;
            double prevCost, newCost;
            prevCost = newCost = slnDict[ rootInd ][ derivNum ].cost;
            MultiOpt::initKBest();

            while ( derivNum < totalDerivs && prevCost == newCost ) {
                MultiOpt::lazyKthBest( H, tinfo, penalty, slnDict, rootInd, derivNum, derivNum );
                cerr << "returned from call to lazy" << derivNum << "best!!\n";
                auto d = slnDict[ rootInd ][ derivNum ];
                cout << "DERIVATION # " << derivNum << ": for [(" << tree->getNodeName(rootKey.u()) << ", " << tree->getNodeName(rootKey.v()) <<  ") : (" << rootKey.f() << rootKey.r() << ")] (" << d.cost << ")\n";

                string fname = "outputDir/output.txt." + lexical_cast<string>(derivNum);
                std::ofstream output( fname, std::ios::out | std::ios::app );
                if ( output.is_open() ) {
                    //output.seekp( 0 );
                    for ( auto f : d.flips ) {
                        cout << " [ " << tree->getNodeName(get<0>(f)) << ", " << tree->getNodeName(get<1>(f)) << " : " << get<2>(f) <<  "]\n";
                        output << tree->getNodeName(get<0>(f)) << "\t" << tree->getNodeName(get<1>(f)) << "\t" << get<2>(f).substr(0,1) << "\n";
                    }
                    output.close();
                } else {
                    throw std::runtime_error("could not open file ("+fname+")");
                }

                prevCost = newCost;
                newCost = d.cost;

                derivNum++;
            }

            delete tree;
        } catch (const Exception &e) {
            cout << "Error when reading tree : " << e.what() << endl;
        }
        delete newickReader;

    } catch(exception& e) {

        cerr << "Caught Exception: [" << e.what() << "]\n";
        abort();

    }

    return 0;
}
