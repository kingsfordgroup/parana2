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
#include <cstdlib>

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
#include "cpplog.hpp"

/** CLN **/
#include <cln/io.h>
#include <cln/rational_io.h>
#include <cln/integer_io.h>


// The logger
#include "cpplog.hpp"
#include "model.hpp"

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

#ifdef CPPLOG_SYSTEM_IDS
#include <boost/interprocess/detail/os_thread_functions.hpp>
using namespace boost::interprocess::detail;
#endif

int main( int argc, char **argv ){
    // initialize the logger
    cpplog::FileLogger log( "log.txt"  );

    po::options_description desc("Allowed Options");

    desc.add_options()
        ("help,h", "produce this message")
        ("target,t", po::value< string >(), "extant graph file")
        ("ratio,r", po::value< double >()->default_value(1.0), "ratio of creation to deletion cost")
        ("beta,b", po::value< double >()->default_value(60.0), "scale factor for cost classes")
        ("undir,u", po::value< bool >()->zero_tokens() , "graph is undirected")
        ("output,o", po::value< string >()->default_value("edgeProbs.txt"), "output file containing edge probabilities")
        ("numOpt,k", po::value< size_t>()->default_value(10), "number of near-optimal score classes to use")
        ("model,m", po::value< string >()->default_value(""), "file containing probabilistic model")
        ("timePenalty,p", po::value< double >()->default_value(0.0), "amount to penalize flips between nodes whose time intervals don't overlap'")
        ("old,x", po::value< bool >()->zero_tokens(), "run using \"old\" algorithm")
        ("prob,p", po::value< bool >()->zero_tokens(), "run using \"max. likelihood\" algorithm")
        ("lazy,l", po::value< bool >()->zero_tokens(), "run using the \"lazy\" algorithm")
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
            for ( const auto& tn : treeNames ) { cerr << "Input tree : " << tn << "\n"; }
        }
        string graphName = ap["target"].as<string>();
        string outputName = ap["output"].as<string>();

        double penalty = ap["timePenalty"].as<double>();
        double creationCost = ap["ratio"].as<double>();
        double deletionCost = 1.0;

        bool undirected = ap.isPresent("undir");
        bool directed = !undirected;
        size_t k = ap["numOpt"].as<size_t>();

        LOG_INFO(log) << "Input graph is " << (undirected  ? "Undirected" : "Directed") << "\n";
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
                LOG_INFO(log) << "rerooted\n";
                // relabel all of the nodes
                tree->resetNodesId();
                LOG_INFO(log) << "relabeled\n";

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
                auto& UG = get<undirectedGraphT>(G);
                //Add all leaf nodes to the graph
                for ( auto nameId : leafIDMap ) {
                    if ( nameId.first.find("LOST") == nameId.first.npos  ) {
                        auto u = add_vertex(UG);
                        UG[u].name = nameId.first;
                        UG[u].idx = nameId.second;
                    }
                }


                if ( isAdjList ) {
                    GraphUtils::readFromAdjacencyList( graphName, leafIDMap, UG );
                } else {
                    GraphUtils::readFromMultilineAdjacencyList(graphName, leafIDMap, UG );
                }

            } else {
                G = directedGraphT();
                auto& DG = get<directedGraphT>(G);
                // Add all leaf nodes to the graph
                for ( auto nameId : leafIDMap ) {
                    if ( nameId.first.find("LOST") == nameId.first.npos  ) {
                        auto u = add_vertex(DG);
                        DG[u].name = nameId.first;
                        DG[u].idx = nameId.second;
                    }
                }

                GraphUtils::readFromMultilineAdjacencyList(graphName, leafIDMap, DG);
            }

            //typedef graph_traits < graphT >::vertex_descriptor vertexDescT;
            //typedef graph_traits < graphT >::edge_descriptor
            //edgeDescT;
            if ( ap.isPresent("prob") ) {

                LOG_INFO(log) << "Max. likelihood algorithm\n";
                string modelFile = ap["model"].as<string>();

                LOG_INFO(log) << "Reading model from file " << modelFile << "\n";
                Model model( modelFile, tinfo, tree );

                auto H = MultiOpt::buildMLSolutionSpaceGraph( tree, tinfo, model, directed );

                vector<size_t> order; order.reserve( H->order() );
                MultiOpt::topologicalOrder( H, order );


                auto rootKey = FlipKey( tree->getRootId(), tree->getRootId(), false, false );
                auto rootInd = H->index(rootKey);

                //vector<size_t> order; order.reserve( H->order() );
                //MultiOpt::topologicalOrderQueue( H, rootInd, order );

                MultiOpt::slnDictT slnDict;
                if ( undirected ) {
                    MultiOpt::MLLeafCostDict( H, tree, get<undirectedGraphT>(G), directed, creationCost, deletionCost, slnDict);
                } else {
                    MultiOpt::MLLeafCostDict( H, tree, get<directedGraphT>(G), directed, creationCost, deletionCost, slnDict);
                }

                MultiOpt::probabilistic<CostClass<EdgeDerivInfoEager>>(H, model, tree, order, slnDict, outputName, keyList );
            } else {

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
            auto beta =  ap["beta"].as<double>();

            if (ap.isPresent("lazy")) {
                // lazy
                auto derivs = MultiOpt::initKBest( H, order, slnDict );
                std::vector<size_t> edges = MultiOpt::viterbiPass(H, derivs, order);

                auto vstr = [&]( const FlipKey& vert ) -> string {
                    auto uname = tree->getNodeName(vert.u());
                    auto vname = tree->getNodeName(vert.v());
                    if (uname > vname) {
                        auto tmp = vname;
                        vname = uname;
                        uname = tmp;
                    }
                    auto fstr = vert.f() ? "true" : "false";
                    auto rstr = vert.r() ? "true" : "false";
                    return uname+"\t"+vname+"\t"+fstr+"\t"+rstr;
                };

                vector<size_t> q;
                q.push_back(rootInd);
                while ( q.size() > 0 ) {
                    auto vit = q.back();
                    auto eid = edges[vit];
                    auto e = H->edge(eid);
                    auto isInternal = H->incident(vit).size() > 0;
                    q.pop_back();
                    cerr << "satisfying " << vstr(H->vertex(vit)) << " using [ ";
                    for ( auto tid : e.tail() ) {
                        cerr << vstr( H->vertex(tid) ) << ", ";
                        if (isInternal) { q.push_back(tid);}
                    }
                    cerr << "]\n";
                }

                exit(0);
                MultiOpt::lazyKthBest( H, rootInd, k, k, derivs );
                MultiOpt::computePosteriors(H, tree, order, derivs, outputName, keyList, beta);
            } else {
                if ( ap.isPresent("old") ) {
                    LOG_INFO(log) << "Old algo\n";
                    MultiOpt::viterbiCount(H, tree, tinfo, penalty, order, slnDict, countDict, k, outputName, keyList, beta);
                } else {
                    LOG_INFO(log) << "New algo\n";
                    // eager
                    MultiOpt::viterbiCountNew<CostClass<EdgeDerivInfoEager>>(H, tree, tinfo, penalty, order, slnDict, countDict, k, outputName, keyList, beta);
                }
            }
        }
        } catch (const Exception &e) {
            LOG_ERROR(log) << "Error when reading tree : " << e.what() << endl;
        }
    } catch(exception& e) {
        LOG_ERROR(log) << "Caught Exception: [" << e.what() << "]\n";
        abort();
    }

    return 0;
}
