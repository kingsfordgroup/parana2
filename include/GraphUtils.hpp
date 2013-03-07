#ifndef GRAPHUTILS_HPP
#define GRAPHUTILS_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/variant.hpp>
#include <boost/lexical_cast.hpp>

namespace GraphUtils {
    using std::cout;
    using std::string;
    using std::ifstream;
    using std::getline;
    using std::vector;
    using std::unordered_map;
    using boost::split;
    using boost::is_any_of;
    using boost::token_compress_on;
    using boost::adjacency_list;
    using boost::graph_traits;
    using boost::vecS;
    using boost::listS;
    using boost::directedS;
    using boost::lexical_cast;
    using std::for_each;
    using boost::variant;
    using boost::add_edge;

    struct Node {
        int idx;
        string name;
    };

    struct Edge {
        double weight;
    };


    typedef adjacency_list<boost::hash_setS, vecS, boost::directedS, GraphUtils::Node, GraphUtils::Edge > directedGraphT;

    typedef adjacency_list<boost::hash_setS, vecS, boost::undirectedS, GraphUtils::Node, GraphUtils::Edge > undirectedGraphT;

    template< typename GraphT >
    struct TreeGraphPair {
        int treeID;
        typename boost::graph_traits< GraphT >::vertex_descriptor graphVert;
    };


    template<typename GraphT>
    void readFromAdjacencyList( const string& fname, GraphT& G ) {
        typedef typename boost::graph_traits< GraphT >::vertex_descriptor Vertex;
        typedef typename boost::graph_traits< GraphT >::edge_descriptor Edge;
        typedef unordered_map<string,Vertex> svMap;

        svMap namePosMap;
        bool inserted;
        Vertex u,v;
        typename unordered_map<string,Vertex>::iterator pos;

        string line;
        typedef vector< string > splitVectorT;
        ifstream gfile(fname);
        size_t numInsertedVerts = 0;
        if ( gfile.is_open() ) {
            while( gfile.good() ) {
                getline( gfile, line, '\n' );
                if ( gfile.eof() ) { break; }
                boost::algorithm::trim(line);
                auto vline = line.substr( 0, line.find_first_of('#') );
                splitVectorT splitVec;
                split( splitVec, vline, is_any_of(" \t"), token_compress_on );

                if ( splitVec.size() > 0  and vline.size() > 0 ) {
                    auto fromVert = splitVec[0];
                    boost::tie( pos, inserted ) = namePosMap.insert( std::make_pair(fromVert,Vertex()) );
                    if (inserted) {
                        ++numInsertedVerts;
                        u = add_vertex(G);
                        G[u].name = fromVert;
                        // This will happen later
                        // G[u].idx = nameMap[fromVert];
                        pos->second = u;
                    } else {
                        u = pos->second;
                    }

                    for( auto tgtIt = splitVec.begin() + 1; tgtIt != splitVec.end(); ++tgtIt ) {
                        auto& tgt = *tgtIt;
                        boost::tie(pos, inserted) = namePosMap.insert(std::make_pair(tgt, Vertex()));
                        if (inserted) {
                            ++numInsertedVerts;
                            v = add_vertex(G);
                            G[v].name = tgt;
                            // This will happen later
                            // G[v].idx = nameMap[tgt];
                            pos->second = v;
                        } else {
                            v = pos->second;
                        }

                        Edge e; bool i;
                        boost::tie(e,i) = add_edge(u,v,G); 
                        G[e].weight = 1.0;
                    } 
                }

            }
            
            if ( namePosMap.size() != boost::num_vertices(G) ) {
                std::cerr << "(namePosMap.size() = " << namePosMap.size() << ") != ("
                          << "(order(G) = " << boost::num_vertices(G) << ") : Error building the graph, aborting\n";
                std::abort();
            }
        }
        gfile.close();

    }

    template<typename GraphT>
    void readFromMultilineAdjacencyList( const string& fname, GraphT& G ) {
        typedef typename boost::graph_traits< GraphT >::vertex_descriptor Vertex;
        typedef typename boost::graph_traits< GraphT >::edge_descriptor Edge;

        typedef unordered_map<string,Vertex> svMap;

        svMap namePosMap;
        bool inserted;
        Vertex u,v;
        typename unordered_map<string,Vertex>::iterator pos;

        bool headLine = false;
        size_t remEdgeLine = 0;
        string line;
        typedef vector< string > splitVectorT;
        ifstream gfile(fname);
        if ( gfile.is_open() ) {
            while( gfile.good() ) {
                getline( gfile, line, '\n' );
                if ( gfile.eof() ) { break; }
                boost::algorithm::trim(line);
                auto vline = line.substr( 0, line.find_first_of('#') );

                if (vline.length() == 0) { continue; }

                splitVectorT splitVec;
                split( splitVec, vline, is_any_of(" \t"), token_compress_on );

                if ( splitVec.size() > 0  and vline.size() > 0 ) {
                    auto fromVert = splitVec[0];
                    remEdgeLine = lexical_cast<size_t>(splitVec[1]);

                    boost::tie( pos, inserted ) = namePosMap.insert( std::make_pair(fromVert,Vertex()) );
                    if (inserted) {
                        u = add_vertex(G);
                        G[u].name = fromVert;
                        // This will happen later
                        // G[u].idx = nameMap[fromVert];
                        pos->second = u;
                    } else {
                        u = pos->second;
                    }

                    while ( remEdgeLine > 0 ) {

                        getline( gfile, line, '\n' );
                        boost::algorithm::trim(line);
                        vline = line.substr( 0, line.find_first_of('#') );
                        split( splitVec, vline, is_any_of(" \t"), token_compress_on );

                        auto toVert = splitVec[0];
                        double weight = lexical_cast<double>(splitVec[1]);


                        boost::tie(pos, inserted) = namePosMap.insert(std::make_pair(toVert, Vertex()));
                        if (inserted) {
                            v = add_vertex(G);
                            G[v].name = toVert;
                            // This will happen later
                            // G[v].idx = nameMap[toVert];
                            pos->second = v;
                        } else {
                            v = pos->second;
                        }

                        Edge e; bool i;
                        boost::tie(e,i) = add_edge(u,v,G);
                        G[e].weight = weight;
                        remEdgeLine--;
                    }
                }

            }

            if ( namePosMap.size() != boost::num_vertices(G) ) {
                std::cerr << "(namePosMap.size() = " << namePosMap.size() << ") != ("
                          << "(order(G) = " << boost::num_vertices(G) << ") : Error building the graph, aborting\n";
                std::abort();
            }
        }
        gfile.close();

    }





}

#endif // GRAPHUTILS_HPP
