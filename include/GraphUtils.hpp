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
#include "utils.hpp"

namespace GraphUtils {
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

    template<typename GraphT>
    void readFromAdjacencyList( const string& fname, unordered_map<string,int>& nameMap, GraphT& G ) {
        typedef typename boost::graph_traits< GraphT >::vertex_descriptor Vertex;
        typedef typename boost::graph_traits< GraphT >::edge_descriptor Edge;
        typedef unordered_map<string,Vertex> svMap;

        for (auto e : nameMap) { cout << "nameMap[" << e.first << "] = " << e.second << "\n";}

        svMap namePosMap;
        bool inserted;
        Vertex u,v;
        typename unordered_map<string,Vertex>::iterator pos;

        string line;
        typedef vector< string > splitVectorT;
        ifstream gfile(fname);
        if ( gfile.is_open() ) {
            while( gfile.good() ) {
                getline( gfile, line, '\n' );
                if ( gfile.eof() ) { break; }
                boost::algorithm::trim(line);
                cout << "LINE = " << line << "\n";
                auto vline = line.substr( 0, line.find_first_of('#') );
                splitVectorT splitVec;
                split( splitVec, vline, is_any_of(" \t"), token_compress_on );
                cout << "TOKENS " << "\n";
                for_each( splitVec.begin(), splitVec.end(), [](const string& s){ cout << " " << s << "\n"; });
                cout << "\n";

                if ( splitVec.size() > 0 ) {
                    auto fromVert = splitVec[0];
                    boost::tie( pos, inserted ) = namePosMap.insert( std::make_pair(fromVert,Vertex()) );
                    if (inserted) {
                        u = add_vertex(G);
                        G[u].name = fromVert;
                        G[u].idx = nameMap[fromVert];
                        pos->second = u;
                    } else {
                        u = pos->second;
                    }
                    //cerr << "G[ " << u << "] :: " << " name = " << G[u].name << ", id = " << G[u].idx << "\n";

                    //cout << "source: " << fromVert << " :: ";
                    auto tgtIt = splitVec.begin(); tgtIt++;
                    for_each( tgtIt, splitVec.end(), [&](const string& tgt) {

                            boost::tie(pos, inserted) = namePosMap.insert(std::make_pair(tgt, Vertex()));
                            if (inserted) {
                                v = add_vertex(G);
                                G[v].name = tgt;
                                G[v].idx = nameMap[tgt];
                                pos->second = v;
                            } else {
                                v = pos->second;
                            }

                            Edge e; bool i;
                            boost::tie(e,i) = add_edge(u,v,G);
                            G[e].weight = 1.0;
                        } );
                    //cout << "\n";
                }

            }
        }
        gfile.close();

    }

    template<typename GraphT>
    void readFromMultilineAdjacencyList( const string& fname, unordered_map<string,int>& nameMap, GraphT& G ) {
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
                for_each( splitVec.begin(), splitVec.end(), [](const string& s){ cout << " " << s << "\n"; });

                if ( splitVec.size() > 0 ) {
                    auto fromVert = splitVec[0];
                    remEdgeLine = lexical_cast<size_t>(splitVec[1]);

                    boost::tie( pos, inserted ) = namePosMap.insert( std::make_pair(fromVert,Vertex()) );
                    if (inserted) {
                        u = add_vertex(G);
                        G[u].name = fromVert;
                        G[u].idx = nameMap[fromVert];
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
                            G[v].idx = nameMap[toVert];
                            pos->second = v;
                        } else {
                            v = pos->second;
                        }

                        Edge e; bool i;
                        boost::tie(e,i) = add_edge(u,v,G);
                        G[e].weight = weight;
                        cerr << "added edge " << G[u].name << ", " << G[v].name << " : " << G[e].weight << "\n";
                        remEdgeLine--;
                    }
                }

            }
        }
        gfile.close();

    }





}
