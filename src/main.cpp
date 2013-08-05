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
#include "CrossValidationParser.hpp"

/** Bio++ **/
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Io/NexusIoTree.h>
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
                 parana2 {pars,prob} [opts]
                 parana2 -h | --help
                 parana2 -v | --version

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
            std::cerr << "\nFor information on a specific method, try \"parana2 {pars,prob} --help\"\n";
            std::exit(1);
        }
    }
}

int main( int argc, char **argv ) {
    // initialize the logger
    cpplog::FileLogger log( "log.txt"  );
    // stderr logger
    cpplog::StdErrLogger slog;

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
    ("cross,c", po::value< string >(), "perform a cross-validation on the extant edges")
    ("target,t", po::value< string >(), "extant graph file")
    ("undir,u", po::value< bool >()->zero_tokens() , "graph is undirected")
    ("output,o", po::value< string >()->default_value("edgeProbs.txt"), "output file containing edge probabilities")
    ("dupHist,d", po::value< vector<string> >(), "duplication history input file(s) [ either 1 or 2 newick format trees ]")
    ;

    // Those options only relevant to the parsimony method
    po::options_description parsimony("Parsimony Specific Options");
    parsimony.add_options()
    ("del", po::value< double >()->default_value(1.0), "deletion cost")
    ("ratio,r", po::value< double >()->default_value(1.2), "ratio of creation to deletion cost")
    ("beta,b", po::value< double >()->default_value(60.0), "scale factor for cost classes")
    ("numOpt,k", po::value< size_t>()->default_value(40), "number of near-optimal score classes to use")
    ("timePenalty,p", po::value< double >()->default_value(1.0), "amount to penalize flips between nodes whose time intervals don't overlap'")
    ("single,s", po::value<bool>()->zero_tokens(), "compute a single optimal set of flips (i.e. \"Parana 1\")")
    ("old,x", po::value< bool >()->zero_tokens(), "run using \"old\" algorithm (deprecated)")
    ("lazy,l", po::value< bool >()->zero_tokens(), "run using the \"lazy\" algorithm (not yet implemented)")
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
            LOG_INFO(slog) << "only 1 or 2 duplication histories can be provided; you provided " + treeNames.size() << "\n";
            abort();
        } else {
            for ( const auto & tn : treeNames ) {
                LOG_INFO(slog) << "Input tree : " << tn << "\n";
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

            vector<size_t> keyList;
            if ( treeNames.size() > 1 ) {

                auto r1 = tree->getSonsId(rId)[0];
                auto r2 = tree->getSonsId(rId)[1];

                keyList.push_back(r1);
                keyList.push_back(r2);
            }
            
           
            tinfo.extantInterval[ rId ] = Utils::ExistenceInterval( -std::numeric_limits<double>::infinity(), 0.0 );
            Utils::Trees::prepareTree( tree, tinfo, rId );

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

            // unordered_map<string,int> leafIDMap;
            // for ( auto nid : tree->getLeavesId() ) {
            //     leafIDMap[ Utils::Trees::getName(tree, nid) ] = nid;
            // }

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

                std::cerr << "reading graph ... ";
                if ( isAdjList ) { // The target graph is unweighted (binary)
                    GraphUtils::readFromAdjacencyList( graphName, UG );
                } else { // The target graph has edge weights
                    GraphUtils::readFromMultilineAdjacencyList(graphName, UG );
                }
                std::cerr << "done\n";

                // Augment the graph nodes with the correpsonding
                // node ids from the phylogeny
                std::cerr << "augmenting graph ... ";
                auto vp = boost::vertices(UG);
                unordered_set<int> leavesCovered;
                for ( auto it = vp.first; it != vp.second; ++it ) {
                    auto lid = tree->getLeafId( UG[*it].name );
                    UG[*it].idx = lid;
                    leavesCovered.insert(lid);
                }
                std::cerr << "done\n";
                if ( leavesCovered.size() != boost::num_vertices(UG) ) {
                    std::cerr << "I added " << leavesCovered.size() << " IDs to the graph, but "
                              << "there were " << boost::num_vertices(UG) << " graph nodes.\n";
                    std::cerr << "These numbers should be the same; aborting\n";
                    std::abort();
                }

            } else { // The target network is directed
                G = directedGraphT();
                auto &DG = get<directedGraphT>(G);

                if ( isAdjList ) { // The target graph is unweighted (binary)
                    GraphUtils::readFromAdjacencyList( graphName, DG );
                } else { // The target graph has edge weights
                    GraphUtils::readFromMultilineAdjacencyList(graphName, DG );
                }

                // Augment the graph nodes with the correpsonding
                // node ids from the phylogeny
                auto vp = boost::vertices(DG);
                unordered_set<int> leavesCovered;
                for ( auto it = vp.first; it != vp.second; ++it ) {
                    auto lid = tree->getLeafId( DG[*it].name );
                    DG[*it].idx = lid;
                    leavesCovered.insert(lid);
                }
                if ( leavesCovered.size() != boost::num_vertices(DG) ) {
                    std::cerr << "I added " << leavesCovered.size() << " IDs to the graph, but "
                              << "there were " << boost::num_vertices(DG) << " graph nodes.\n";
                    std::cerr << "These numbers should be the same; aborting\n";
                    std::abort();
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
                    MultiOpt::MLLeafCostDict( H, tree, tinfo, get<undirectedGraphT>(G), directed, creationCost, deletionCost, slnDict);
                } else {
                    MultiOpt::MLLeafCostDict( H, tree, tinfo, get<directedGraphT>(G), directed, creationCost, deletionCost, slnDict);
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


                auto performReconstruction = [&]( const std::string& outputFileName, MultiOpt::slnDictT& slnDict ) -> bool {
                    //auto rootKey = FlipKey( tree->getRootId(), tree->getRootId(), false, false, true, true );
                    //auto rootInd = H->index(rootKey);

                    // Count the # of opt slns.
                    MultiOpt::countDictT countDict;
                    //countDict.set_empty_key(-1);
                    countDict.resize(order.size());
                
                    auto beta =  ap["beta"].as<double>();

                    if (ap.isPresent("single")) {
                        MultiOpt::viterbi<CostClass<EdgeDerivInfoEager>>(H, tree, tinfo, penalty, order, slnDict, outputFileName);
                    } else if (ap.isPresent("lazy")) {
                        // lazy
                        auto derivs = MultiOpt::initKBest( H, order, slnDict );
                        std::vector<size_t> edges = MultiOpt::viterbiPass(H, derivs, order);

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
                            std::cerr << "KEY LIST SIZE IS " << keyList.size() << "\n";
                            // eager
                            MultiOpt::viterbiCountNew<CostClass<EdgeDerivInfoEager>>(H, tree, tinfo, penalty, order, slnDict, countDict, k, outputFileName, keyList, beta);
                            std::cerr << "DONE" << std::endl;
                            LOG_INFO(log) << "DONE\n";
                        }
                    }
                }; // end of reconstruction function

            auto doCrossValidation = ap.isPresent("cross");
            string crossValidationFile;
            if ( doCrossValidation ) {
                crossValidationFile = ap["cross"].as<string>();
            }
            path outputPath(outputName);
            switch (doCrossValidation) {
                case true:
                  {
                  CrossValidationTestParser parser( crossValidationFile );
                  auto crossValTests = parser.parse();

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

                    auto isLostNode = [&]( const std::string& n ) -> bool { return n.find("LOST") != n.npos; };
                    //
                    //  We're only interested in potential "extant" edges, so only output those
                    //
                    unordered_set<size_t> lostNodes;
                    unordered_set<string> vertNames;
                    auto vp = boost::vertices(Graph);
                    for ( auto it = vp.first; it != vp.second; ++it ) {
                        if (isLostNode( Graph[*it].name )) { lostNodes.insert(Graph[*it].idx); }
                        vertNames.insert(Graph[*it].name);
                    }

                    for ( auto it = vp.first; it != vp.second; ++it ) {
                        auto u = *it;
                        auto idxU = Graph[u].idx;
                        if( isLostNode(Graph[u].name) ) { continue; }

                        for ( auto it2 = it; it2 != vp.second; ++it2 ) {
                            auto v = *it2;
                            auto idxV = Graph[v].idx;
                            if( isLostNode(Graph[v].name) ) { continue; }

                            auto i1 = std::min(idxU, idxV);
                            auto i2 = std::max(idxU, idxV);
                            
                            if ( Utils::Trees::sameSpecies( tree, i1, i2 ) ) {
                                keyList.emplace_back( H->index(FlipKey(i1, i2, FlipState::both)) );
                            }
                        }
                    }

                    // Build a name -> vertex map for easy access of edges
                    unordered_map<string, VertexT> nameVertMap;
                    for ( auto it = vp.first; it != vp.second; ++it ) {
                        nameVertMap[Graph[*it].name] = *it;
                    }

     
                    // LEAF COST DICT
                    MultiOpt::slnDictT slnDict;
                    slnDict.resize(order.size());
                    if ( undirected ) {
                        MultiOpt::leafCostDict( H, tree, tinfo, get<undirectedGraphT>(G), directed, creationCost, deletionCost, slnDict);
                    } else {
                        MultiOpt::leafCostDict( H, tree, tinfo, get<directedGraphT>(G), directed, creationCost, deletionCost, slnDict);
                    }
                    // DONE LEAF COST DICT

                    for ( auto& cvtest : crossValTests.cvSets ) {
                        
                        std::vector<double> weights;
                        weights.reserve(boost::num_edges(Graph));

                        // ------- Remove the edges from the graph
                        for ( auto& stredge : cvtest.edges ) {
                            auto u = nameVertMap[stredge.u];
                            auto v = nameVertMap[stredge.v];

                            EdgeT edge; bool exists;
                            std::tie(edge,exists) = boost::edge(u,v,Graph);
                            weights.emplace_back(Graph[edge].weight);
                            boost::remove_edge(u,v,Graph);

                            auto idxU = Graph[u].idx; auto idxV = Graph[v].idx;
                            if ( idxU > idxV ) { std::swap(idxU, idxV); }
                            auto edgeKey = H->index(FlipKey( idxU, idxV, FlipState::both ));
                            auto noEdgeKey = H->index(FlipKey( idxU, idxV, FlipState::none ));
                            // it no longer exists so this costs something
                            auto prevExistCost = slnDict[edgeKey][0].cost;
                            auto prevNoExistCost = slnDict[noEdgeKey][0].cost;
                            if ( prevExistCost > 0.0 or prevNoExistCost == 0.0 ) {
                                std::cerr << "WTF?!\n";
                                std::abort();
                            }
                            slnDict[edgeKey][0].cost = deletionCost;
                            // this is now free
                            slnDict[noEdgeKey][0].cost = 0.0;
                        }
                        std::cerr << "\n";
                        // ------- Done removing edges
                        
                        string reconFileName = boost::str( boost::format("%1%/%2%.txt") % outputName % cvtest.name );
                        performReconstruction( reconFileName, slnDict );

                        // ------- Put the edges back into the graph
                        size_t weightIdx = 0;
                        for ( auto& stredge : cvtest.edges ) {
                            auto u = nameVertMap[stredge.u];
                            auto v = nameVertMap[stredge.v];
                            auto weight = weights[weightIdx];
                            ++weightIdx;
                            // Add the edge back to the graph
                            EdgeT edge; bool exists;
                            std::tie(edge, exists) = boost::add_edge(u,v,Graph);
                            Graph[edge].weight = weight;

                            auto idxU = Graph[u].idx; auto idxV = Graph[v].idx;
                            if ( idxU > idxV ) { std::swap(idxU, idxV); }
                            auto edgeKey = H->index(FlipKey( idxU, idxV, FlipState::both ));
                            auto noEdgeKey = H->index(FlipKey( idxU, idxV, FlipState::none ));
                            // Put the leaf cost dictionary back to the original state
                            slnDict[edgeKey][0].cost = 0.0;
                            slnDict[noEdgeKey][0].cost = creationCost;

                        }   
                        // ------- Done putting back
                    }

                    // ************ ORIGINAL CROSS VALIDATION CODE ************ 
                    // VertexT u,v;
                    // for( auto& uv : edges ) {
                    //     std::tie(u,v) = uv;
                    //     EdgeT edge; bool exists;
                    //     std::tie(edge,exists) = boost::edge(u,v,Graph);
                    //     auto weight = Graph[edge].weight;
                    //     // Remove the edge from the graph
                    //     boost::remove_edge(u,v,Graph);

                    //     // Modify the leaf cost dictionary
                    //     auto idxU = Graph[u].idx;
                    //     auto idxV = Graph[v].idx;
                    //     if ( idxU > idxV ) { std::swap(idxU, idxV); }
                    //     auto edgeKey = H->index(FlipKey( idxU, idxV, FlipState::both ));
                    //     auto noEdgeKey = H->index(FlipKey( idxU, idxV, FlipState::none ));

                    //     // it no longer exists so this costs something
                    //     slnDict[edgeKey][0].cost = deletionCost;
                    //     // this is now free
                    //     slnDict[noEdgeKey][0].cost = 0.0;

                    //     // Perform the reconstruction
                    //     auto nameU = Graph[u].name; auto nameV = Graph[v].name;
                    //     auto reconFileName = boost::str( boost::format("%1%/removed@%2%#%3%@txt") % outputName % nameU % nameV );
                    //     performReconstruction( reconFileName, slnDict );

                    //     // Add the edge back to the graph
                    //     std::tie(edge, exists) = boost::add_edge(u,v,Graph);
                    //     Graph[edge].weight = weight;

                    //     // Put the leaf cost dictionary back to the original state
                    //     slnDict[edgeKey][0].cost = 0.0;
                    //     slnDict[noEdgeKey][0].cost = creationCost;

                    // }
                    // ************ END OF ORIGINAL CROSS VALIDATION CODE ************ 

                  } else { // The input graph is directed
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

                    // LEAF COST DICT
                    MultiOpt::slnDictT slnDict;
                    slnDict.resize(order.size());
                    if ( undirected ) {
                        MultiOpt::leafCostDict( H, tree, tinfo, get<undirectedGraphT>(G), directed, creationCost, deletionCost, slnDict);
                    } else {
                        MultiOpt::leafCostDict( H, tree, tinfo, get<directedGraphT>(G), directed, creationCost, deletionCost, slnDict);
                    }
                    // DONE LEAF COST DICT

                    /** original Leave-one-out cross-validation **/
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
                        performReconstruction( reconFileName, slnDict );

                        // Add the edge back to the graph
                        std::tie(edge, exists) = boost::add_edge(u,v,Graph);
                        Graph[edge].weight = weight;
                    }


                  }
                  }
                  break;

                case false:
                    MultiOpt::slnDictT slnDict;
                    slnDict.resize(order.size());
                    if ( undirected ) {
                        MultiOpt::leafCostDict( H, tree, tinfo, get<undirectedGraphT>(G), directed, creationCost, deletionCost, slnDict);
                    } else {
                        MultiOpt::leafCostDict( H, tree, tinfo, get<directedGraphT>(G), directed, creationCost, deletionCost, slnDict);
                    }

                    performReconstruction( outputName, slnDict );
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
