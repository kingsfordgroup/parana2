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
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/progress.hpp>
#include <boost/format.hpp>
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
#include <algorithm>
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
#include "PhyloXMLParser.hpp"

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
using boost::filesystem::path;
using boost::filesystem::exists;
using boost::filesystem::is_directory;
using boost::filesystem::create_directories;
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

void printUsage() {
    auto usage = R"(Usage:
                 parana++ {pars,prob} [opts]
                 parana++ -h | --help
                 parana++ -v | --version

                 Arguments:
                 pars    Use the parsimonious recovery algorithm
                 prob    Use a probabilistic (ML) recover algorithm
                 )";
    std::cout << usage;
}

void printHelpUsage( ArgParser &ap,
                     po::options_description &general,
                     po::options_description &prob,
                     po::options_description &parsimony) {

    if ( ap.isPresent("help") ) {
        if ( ap.isPresent("method") ) {
            auto method = ap["method"].as<std::string>();
            if ( method == "prob" ) {
                std::cout << general;
                std::cout << prob;
                std::exit(1);
            } else if ( method == "pars" ) {
                std::cout << general;
                std::cout << parsimony;
                std::exit(1);
            } else {
                printUsage();
                std::exit(1);
            }
        } else {
            printUsage();
            std::cerr << "\nFor information on a specific method, try \"parana++ {pars,prob} --help\"\n";
            std::exit(1);
        }
    }
}

int main( int argc, char **argv ) {
    // initialize the logger
    cpplog::FileLogger log( "log.txt"  );

    // The options to choose the method (parsimony | probabilistic)
    po::positional_options_description p;
    p.add("method", 1);
    po::options_description pos("");
    pos.add_options()
    ("method", po::value< string >(), "method to use; one of {prob,pars}");

    // The help option
    po::options_description help("Help Options");
    help.add_options()
    ("help,h", po::value<bool>()->zero_tokens(), "display help for a specific method");

    // The general options are relevant to either method
    po::options_description general("Genearal Options");
    general.add_options()
    ("cross,c", po::value< bool >()->zero_tokens(), "perform a cross-validation on the extant edges")
    ("target,t", po::value< string >(), "extant graph file")
    ("undir,u", po::value< bool >()->zero_tokens() , "graph is undirected")
    ("output,o", po::value< string >()->default_value("edgeProbs.txt"), "output file containing edge probabilities")
    ("dupHist,d", po::value< vector<string> >(), "duplication history input file(s) [ either 1 or 2 newick format trees ]")
    ;

    // Those options only relevant to the parsimony method
    po::options_description parsimony("Parsimony Specific Options");
    parsimony.add_options()
    ("del", po::value< double >()->default_value(1.0), "deletion cost")
    ("ratio,r", po::value< double >()->default_value(1.0), "ratio of creation to deletion cost")
    ("beta,b", po::value< double >()->default_value(60.0), "scale factor for cost classes")
    ("numOpt,k", po::value< size_t>()->default_value(10), "number of near-optimal score classes to use")
    ("timePenalty,p", po::value< double >()->default_value(0.0), "amount to penalize flips between nodes whose time intervals don't overlap'")
    ("old,x", po::value< bool >()->zero_tokens(), "run using \"old\" algorithm")
    ("lazy,l", po::value< bool >()->zero_tokens(), "run using the \"lazy\" algorithm")
    ;

    // Those options only relevant to the probabilistic method
    po::options_description prob("Probabilistic Specific Options");
    prob.add_options()
    ("model,m", po::value< string >()->default_value(""), "file containing probabilistic model")
    ;

    // The options description to parse all of the options
    po::options_description all("Allowed options");
    all.add(help).add(general).add(pos).add(prob).add(parsimony);

    try {
        ArgParser ap( all, p );
        // If there are no extra command line arguments, then print the useage and exit
        if ( argc == 1 ) {
            printUsage();
            std::exit(1);
        }
        ap.parse_args( argc, argv );

        // For any combination of cmd line arguments that won't actually run the program
        // print the appropriate message and exit
        printHelpUsage(ap, general, prob, parsimony);
        auto method = ap["method"].as<std::string>();

        // Get the duplication histories
        vector<string> treeNames = ap["dupHist"].as< vector<string> >();
        // We can't have more than 2
        if ( treeNames.size() > 2 ) {
            cerr << "only 1 or 2 duplication histories can be provided; you provided " + treeNames.size() << "\n";
            abort();
        } else {
            for ( const auto & tn : treeNames ) {
                cerr << "Input tree : " << tn << "\n";
            }
        }

        string graphName = ap["target"].as<string>();
        string outputName = ap["output"].as<string>();

        bool undirected = ap.isPresent("undir");
        bool directed = !undirected;

        double penalty = ap["timePenalty"].as<double>();
        double deletionCost = ap["del"].as<double>();
        double creationCost = ap["ratio"].as<double>() * deletionCost;
        size_t k = ap["numOpt"].as<size_t>();

        std::cerr << "[creation] : [deletion] ratio is [" << creationCost << "] : [" << deletionCost << "]\n";
        std::cerr << "TimePenalty is [" << penalty << "]\n";
        std::cerr << "Considering top " << k << " cost classes\n";

        LOG_INFO(log) << "Input graph is " << (undirected  ? "Undirected" : "Directed") << "\n";
        try {
            PhyloXMLParser parser(treeNames[0]);

            //TreePtrT tree ( Utils::Trees::readNewickTree(treeNames[0]) );
            TreePtrT tree ( parser.parse() );

            
            TreePtrT tree2;

            if ( treeNames.size() > 1 ) {

                PhyloXMLParser t2parser(treeNames[1]);
                auto del = [](TreePtrT::element_type* e) { /*do nothing deleter*/ };
                tree2.reset( t2parser.parse(), del );
                
                auto nodeName = "#preroot#";
                auto node = new bpp::Node(nodeName);

                bpp::BppString specName("#prspec#");
                node->setNodeProperty("S", specName);
                assert( dynamic_cast<bpp::BppString*>(node->getNodeProperty("S"))->toSTL() == "#prspec#" );

                // This is a speciation event
                node->setNodeProperty("D", bpp::BppString("F") );
                assert( dynamic_cast<bpp::BppString*>(node->getNodeProperty("D"))->toSTL() == "F" );

                auto nodeString = new bpp::BppString(nodeName);
                node->setNodeProperty("GN", bpp::BppString(nodeName));
                assert( dynamic_cast<bpp::BppString*>(node->getNodeProperty("GN"))->toSTL() == nodeName );

                auto r1 = tree->getRootNode();
                auto r2 = tree2->getRootNode();

                // relationship with tree 1
                node->addSon(0, tree->getRootNode());
                r1->setDistanceToFather(0.0);

                // relationship with tree 2
                node->addSon(1, tree2->getRootNode());
                r2->setDistanceToFather(0.0);

                // root the tree at our new node
                tree->setRootNode(node);
                std::cerr << "rerooted\n";
                // relabel all of the nodes
                std::cerr << "before relabel tree has " << tree->getNodesId().size() << " nodes\n";
                tree->resetNodesId();
                std::cerr << "relabeled\n";
                std::cerr << "after relabel tree has " << tree->getNodesId().size() << " nodes\n";

                std::cerr << "before\n";
                auto rootId = tree->getRootId();
                std::cerr << "after\n";
                // set the name of the root
                tree->setNodeName( rootId, "#preroot#");
                tree->setVoidBranchLengths(0.0);
                std::cerr << "tree = " << tree << "\n";
                std::cerr << "right before scope exit tree has " << tree->getNodesId().size() << " nodes\n";
            }
            std::cerr << "tree = " << tree << "\n";
            std::cerr << "right after scope exit tree has " << tree->getNodesId().size() << " nodes\n";
            std::cerr << "before label\n";
            std::cerr << "there are " << tree->getNodesId().size() << " nodes\n";
            // put the node names on all the nodes
            Utils::Trees::labelTree(tree);
            std::cerr << "after label\n";

            TreeInfo tinfo("TreeInfo", tree);
            auto rId = tree->getRootId();

            vector<FlipKey> keyList;
            if ( treeNames.size() > 1 ) {

                auto r1 = tree->getSonsId(rId)[0];
                auto r2 = tree->getSonsId(rId)[1];

                keyList.push_back( FlipKey(r1, r2, true, true) );
            }

            tinfo.extantInterval[ rId ] = Utils::ExistenceInterval( -std::numeric_limits<double>::infinity(), 0.0 );
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
            for ( auto nid : tree->getLeavesId() ) {
                leafIDMap[ Utils::Trees::getName(tree, nid) ] = nid;
            }

            // Check if the graph is an adjacency list (this means it's unweighted)
            auto isAdjList = [&]() {
                auto spos = graphName.find_last_of(".");
                auto suffix(graphName.substr( spos ));
                return suffix == ".adj";
            }();

            variant< directedGraphT, undirectedGraphT > G;
            if ( undirected ) { // The target network is undirected
                G = undirectedGraphT();
                auto &UG = get<undirectedGraphT>(G);

                //Add all leaf nodes to the graph
                for ( auto nameId : leafIDMap ) {
                    if ( nameId.first.find("LOST") == nameId.first.npos  ) {
                        auto u = add_vertex(UG);
                        UG[u].name = nameId.first;
                        UG[u].idx = nameId.second;
                    }
                }

                if ( isAdjList ) { // The target graph is unweighted (binary)
                    GraphUtils::readFromAdjacencyList( graphName, leafIDMap, UG );
                } else { // The target graph has edge weights
                    GraphUtils::readFromMultilineAdjacencyList(graphName, leafIDMap, UG );
                }

            } else { // The target network is directed
                G = directedGraphT();
                auto &DG = get<directedGraphT>(G);

                // Add all leaf nodes to the graph
                for ( auto nameId : leafIDMap ) {
                    if ( nameId.first.find("LOST") == nameId.first.npos  ) {
                        auto u = add_vertex(DG);
                        DG[u].name = nameId.first;
                        DG[u].idx = nameId.second;
                    }
                }

                if ( isAdjList ) { // The target graph is unweighted (binary)
                    GraphUtils::readFromAdjacencyList( graphName, leafIDMap, DG );
                } else { // The target graph has edge weights
                    GraphUtils::readFromMultilineAdjacencyList(graphName, leafIDMap, DG );
                }

            }

            // Using probabilistic (Maximum Likelihood) algo.
            if ( method == "prob" ) {

                LOG_INFO(log) << "Max. likelihood algorithm\n";
                string modelFile = ap["model"].as<string>();

                LOG_INFO(log) << "Reading model from file " << modelFile << "\n";
                Model model( modelFile, tinfo, tree );

                auto H = MultiOpt::buildMLSolutionSpaceGraph( tree, tinfo, model, directed );
                //auto H = MultiOpt::buildSolutionSpaceGraph( tree, tinfo, creationCost, deletionCost, penalty, directed );
                
                auto rootKey = FlipKey( tree->getRootId(), tree->getRootId(), false, false, true, true );
                auto rootInd = H->index(rootKey);

                vector<size_t> order; order.reserve( H->order() );
                MultiOpt::topologicalOrder( H, tree, tinfo, order );

                //vector<size_t> order; order.reserve( H->order() );
                //MultiOpt::topologicalOrderQueue( H, rootInd, order );

                MultiOpt::slnDictT slnDict;
                slnDict.resize(order.size());
                if ( undirected ) {
                    MultiOpt::MLLeafCostDict( H, tree, get<undirectedGraphT>(G), directed, creationCost, deletionCost, slnDict);
                } else {
                    MultiOpt::MLLeafCostDict( H, tree, get<directedGraphT>(G), directed, creationCost, deletionCost, slnDict);
                }

                MultiOpt::probabilistic<CostClass<EdgeDerivInfoEager>>(H, model, tree, order, slnDict, outputName, keyList );

            } else if ( method == "pars" ) { // Using parsimony algo.

                auto H = MultiOpt::buildSolutionSpaceGraph( tree, tinfo, creationCost, 
                                                            deletionCost, penalty, 
                                                            directed, MultiOpt::DerivationType::AllHistories );

                auto rootKey = FlipKey( tree->getRootId(), tree->getRootId(), false, false, true, true );
                auto rootInd = H->index(rootKey);

                vector<size_t> order; order.reserve( H->order() );
                MultiOpt::topologicalOrder( H, tree, tinfo, order );

                auto performReconstruction = [&]( const std::string& outputFileName ) -> bool {
                MultiOpt::slnDictT slnDict;
                slnDict.resize(order.size());

                //slnDict.set_empty_key( std::numeric_limits<size_t>::max() );
                if ( undirected ) {
                    MultiOpt::leafCostDict( H, tree, tinfo, get<undirectedGraphT>(G), directed, creationCost, deletionCost, slnDict);
                } else {
                    MultiOpt::leafCostDict( H, tree, tinfo, get<directedGraphT>(G), directed, creationCost, deletionCost, slnDict);
                }

                //auto rootKey = FlipKey( tree->getRootId(), tree->getRootId(), false, false, true, true );
                //auto rootInd = H->index(rootKey);

                // Count the # of opt slns.
                MultiOpt::countDictT countDict;
                //countDict.set_empty_key(-1);
                countDict.resize(order.size());
                
                auto beta =  ap["beta"].as<double>();

                if (ap.isPresent("lazy")) {
                    // lazy
                    auto derivs = MultiOpt::initKBest( H, order, slnDict );
                    std::vector<size_t> edges = MultiOpt::viterbiPass(H, derivs, order);

                    auto vstr = [&]( const FlipKey & vert ) -> string {
                        auto uname = tree->getNodeName(vert.u());
                        auto vname = tree->getNodeName(vert.v());
                        if (uname > vname) {
                            auto tmp = vname;
                            vname = uname;
                            uname = tmp;
                        }
                        auto fstr = vert.f() ? "true" : "false";
                        auto rstr = vert.r() ? "true" : "false";
                        return uname + "\t" + vname + "\t" + fstr + "\t" + rstr;
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
                            if (isInternal) {
                                q.push_back(tid);
                            }
                        }
                        cerr << "]\n";
                    }

                    exit(0);
                    MultiOpt::lazyKthBest( H, rootInd, k, k, derivs );
                    MultiOpt::computePosteriors(H, tree, order, derivs, outputFileName, keyList, beta);
                } else {
                    if ( ap.isPresent("old") ) {
                        LOG_INFO(log) << "Old algo\n";
                        MultiOpt::viterbiCount(H, tree, tinfo, penalty, order, slnDict, countDict, k, outputFileName, keyList, beta);
                    } else {
                        LOG_INFO(log) << "New algo\n";
                        // eager
                        MultiOpt::viterbiCountNew<CostClass<EdgeDerivInfoEager>>(H, tree, tinfo, penalty, order, slnDict, countDict, k, outputFileName, keyList, beta);
                        std::cerr << "DONE" << std::endl;
                        LOG_INFO(log) << "DONE\n";
                    }
                }
            }; // end of reconstruction function

            auto doCrossValidation = ap.isPresent("cross");
            path outputPath(outputName);
            switch (doCrossValidation) {
                case true:
                  // Check if the cross validation directory exists
                  if ( exists(outputPath) ) {
                    if ( !is_directory(outputPath) ) { 
                        std::cerr << "The path given for cross validation output [" << outputPath << 
                                  "] is not a directory\n";
                    }
                  } else {
                    if (!create_directories(outputPath)) {
                        std::cerr << "Could not create the cross validation output directory [ " <<
                                      outputPath << "]\n";
                    }
                  }

                  if ( undirected ) {
                    auto& Graph = get<undirectedGraphT>(G);
                    typedef typename boost::graph_traits< undirectedGraphT >::edge_descriptor EdgeT;
                    typedef typename boost::graph_traits< undirectedGraphT >::vertex_descriptor VertexT;
                    typedef typename boost::graph_traits< undirectedGraphT >::edge_iterator EdgeItT;
                    EdgeItT eit, eend;
                    std::tie(eit, eend) = boost::edges( Graph );
                    std::vector< std::tuple<VertexT,VertexT> > edges;
                    edges.reserve( std::distance(eit, eend) );
                    for ( eit; eit != eend; ++eit ) {
                        edges.push_back(make_tuple( boost::source(*eit,Graph), boost::target(*eit, Graph) ));
                    }

                    VertexT u,v;
                    for( auto& uv : edges ) {
                        std::tie(u,v) = uv;
                        EdgeT edge; bool exists;
                        std::tie(edge,exists) = boost::edge(u,v,Graph);
                        auto weight = Graph[edge].weight;
                        // Remove the edge from the graph
                        boost::remove_edge(u,v,Graph);

                        // Perform the reconstruction
                        auto nameU = Graph[u].name; auto nameV = Graph[v].name;
                        auto reconFileName = boost::str( boost::format("%1%/removed@%2%#%3%@txt") % outputName % nameU % nameV );
                        performReconstruction( reconFileName );

                        // Add the edge back to the graph
                        std::tie(edge, exists) = boost::add_edge(u,v,Graph);
                        Graph[edge].weight = weight;
                    }

                  } else { // The input graph is undirected
                    auto& Graph = get<directedGraphT>(G);
                    typedef typename boost::graph_traits< directedGraphT >::edge_descriptor EdgeT;
                    typedef typename boost::graph_traits< directedGraphT >::vertex_descriptor VertexT;
                    typedef typename boost::graph_traits< directedGraphT >::edge_iterator EdgeItT;
                    EdgeItT eit, eend;
                    std::tie(eit, eend) = boost::edges( Graph );
                    std::vector< std::tuple<VertexT,VertexT> > edges;
                    edges.reserve( std::distance(eit, eend) );
                    for ( eit; eit != eend; ++eit ) {
                        edges.push_back(make_tuple( boost::source(*eit,Graph), boost::target(*eit, Graph) ));
                    }

                    VertexT u,v;
                    for( auto& uv : edges ) {
                        std::tie(u,v) = uv;
                        EdgeT edge; bool exists;
                        std::tie(edge,exists) = boost::edge(u,v,Graph);
                        auto weight = Graph[edge].weight;
                        // Remove the edge from the graph
                        boost::remove_edge(u,v,Graph);

                        // Perform the reconstruction
                        auto nameU = Graph[u].name; auto nameV = Graph[v].name;
                        auto reconFileName = boost::str( boost::format("%1%/removed@%2%#%3%@txt") % outputName % nameU % nameV );
                        performReconstruction( reconFileName );

                        // Add the edge back to the graph
                        std::tie(edge, exists) = boost::add_edge(u,v,Graph);
                        Graph[edge].weight = weight;
                    }
                  }
                  break;

                case false:
                  performReconstruction( outputName );
                  std::cerr << "HERE\n";
                  std::exit(1);
                  break;
            }
        }

        } catch (const Exception &e) {
            LOG_ERROR(log) << "Error when reading tree : " << e.what() << endl;
        }
    } catch (exception &e) {
        LOG_ERROR(log) << "Caught Exception: [" << e.what() << "]\n";
        abort();
    }

    return 0;
}
