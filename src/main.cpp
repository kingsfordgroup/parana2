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
#include "MultiOpt.hpp"

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
using Utils::Trees::TreePtrT;
using GraphUtils::undirectedGraphT;
using GraphUtils::directedGraphT;

#include <cstdlib>

int main( int argc, char **argv ){
    po::options_description desc("Allowed Options");

    desc.add_options()
        ("help,h", "produce this message")
        ("target,t", po::value< string >(), "extant graph file")
        ("ratio,r", po::value< double >()->default_value(1.0), "ratio of creation to deletion cost")
        ("undir,u", po::value< bool >()->zero_tokens() , "graph is undirected")
        ("output,o", po::value< string >()->default_value("edgeProbs.txt"), "output file containing edge probabilities")
        ("numOpt,k", po::value< size_t>()->default_value(10), "number of near-optimal score classes to use")
        ("timePenalty,p", po::value< double >()->default_value(0.0), "amount to penalize flips between nodes whose time intervals don't overlap'")
        ("old,l", po::value< bool >()->zero_tokens(), "run using \"old\" algorithm")
        ("dupHist,d", po::value< vector<string> >(), "duplication history input file(s) [ either 1 or 2 newick format trees ]")
        ;

    try {
        ArgParser ap( desc );
        ap.parse_args( argc, argv );

        vector<string> treeNames = ap["dupHist"].as< vector<string> >();
        if ( treeNames.size() > 2 ) {
            cerr << "only 1 or 2 duplication histories can be provided; you provided "+treeNames.size() << "\n";
            abort();
        } else {
            for ( const auto& tn : treeNames ) { cout << "input tree : " << tn << "\n"; }
        }
        string graphName = ap["target"].as<string>();
        string outputName = ap["output"].as<string>();

        double penalty = ap["timePenalty"].as<double>();
        double creationCost = ap["ratio"].as<double>();
        double deletionCost = 1.0;

        bool undirected = ap.isPresent("undir");
        bool directed = !undirected;
        size_t k = ap["numOpt"].as<size_t>();

        cout << "UNDIRECTED = " << undirected << "\n";
        try {


            TreePtrT tree ( Utils::Trees::readNewickTree(treeNames[0]) );

            if ( treeNames.size() > 1 ) {
                auto tree2 = Utils::Trees::readNewickTree(treeNames[1]).release();
                auto node = new bpp::Node("#preroot#");
                auto r1 = tree->getRootNode();
                auto r2 = tree2->getRootNode();

                // relationship with tree 1
                node->addSon(0, tree->getRootNode());
                // relationship with tree 2
                node->addSon(1, tree2->getRootNode());

                // root the tree at our new node
                tree->setRootNode(node);
                cout << "rerooted\n";
                // relabel all of the nodes
                tree->resetNodesId();
                cout << "relabeled\n";

                auto rootId = tree->getRootId();
                // set the name of the root
                tree->setNodeName( rootId, "#preroot#");
                tree->setVoidBranchLengths(0.0);
            }

            // put the node names on all the nodes
            Utils::Trees::labelTree(tree);


            TreeInfo tinfo("TreeInfo");
            auto rId = tree->getRootId();

            vector<FlipKey> keyList;
            if ( treeNames.size() > 1 ) {

                auto r1 = tree->getSonsId(rId)[0];
                auto r2 = tree->getSonsId(rId)[1];

                keyList.push_back( FlipKey(r1,r2,true,true) );
            }

            tinfo.extantInterval[ rId ] = make_tuple( -std::numeric_limits<double>::infinity(), 0.0 );
            Utils::Trees::prepareTree( tree, tinfo, rId );
            /*
            auto nodes = vector<int>(tree->getNodesId());
            std::sort( nodes.begin(), nodes.end(),
                [&]( const int& n1, const int& n2 ) -> bool {
                    return tree->getNodeName(n1) < tree->getNodeName(n2);
                }
            );
            for ( auto n : nodes ) {
                cerr << tree->getNodeName(n) << ",";
                for ( auto sn : tinfo.leaves[n] ) {
                    cerr << " " << tree->getNodeName(sn);
                }
                cerr << ",";
                for ( auto en : tinfo.enets[n] ) {
                    cerr << " " << en;
                }
                cerr << "\n";
            }
            exit(0);
            */
            // print out the existence intervals for all of the nodes
            /*
            for ( auto nid : tree->getNodesId() ){
                auto nodeName = Utils::Trees::getName(tree,nid);
                cout << "Node [" << nodeName << "] : existance interval = [" <<
                    get<0>(tinfo.extantInterval[nid]) << ", " << get<1>(tinfo.extantInterval[nid]) << "] \n";
            }
            */

            unordered_map<string, int> leafIDMap;
            for( auto nid : tree->getLeavesId() ) { leafIDMap[ Utils::Trees::getName(tree,nid) ] = nid; }

            auto isAdjList = [&](){
                auto spos = graphName.find_last_of(".");
                auto suffix(graphName.substr( spos ));
                return suffix == ".adj";
            }();

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
            if ( ap.isPresent("old") ) {
                cerr << "Old algo\n";
                MultiOpt::viterbiCount(H, tree, tinfo, penalty, order, slnDict, countDict, k, outputName, keyList);
            } else {
                cerr << "New algo\n";
                MultiOpt::viterbiCountNew(H, tree, tinfo, penalty, order, slnDict, countDict, k, outputName, keyList);
            }

            //cout << "The optimal cost solutions have a cost of " << slnDict[rootInd][0].cost << "\n";
            //cout << "There are " << get<1>(countDict[rootInd][0]) << " optimal solutions ";
            //cout << "Near optimal -- (cost, count) : \n";
            //for ( const auto& der : countDict[rootInd] ) {
            //    cout << " (" << get<0>(der) << ",  " << get<1>(der) << ")\n";
            //}
            //cout << "\n";

            /*
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
            */
        } catch (const Exception &e) {
            cout << "Error when reading tree : " << e.what() << endl;
        }
    } catch(exception& e) {

        cerr << "Caught Exception: [" << e.what() << "]\n";
        abort();

    }

    return 0;
}
