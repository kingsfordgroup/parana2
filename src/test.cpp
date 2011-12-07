#include <unordered_map>
#include <unordered_set>
#include <string>
#include <tuple>
#include <iostream>
#include <array>
#include <vector>
//#include "MultiOpt.hpp"
#include <boost/bimap.hpp>
#include <boost/bimap/vector_of.hpp>
#include <boost/bimap/unordered_set_of.hpp>
#include <functional>
#include "FlipKey.hpp"
#include "Derivation.hpp"
#include <chrono>
#include <random>
#include <iterator>

typedef std::mt19937 MyRNG;  // the Mersenne Twister with a popular choice of parameters
uint32_t seed_val;           // populate somehow

MyRNG rng;                   // e.g. keep one global instance (per thread)



#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/labeled_graph.hpp>
#include <google/dense_hash_set>
#include <google/dense_hash_map>

using google::dense_hash_set;
using google::dense_hash_map;
using std::cerr;
using std::tuple;
using std::make_tuple;
using std::get;
using std::unordered_map;
using std::unordered_set;
using std::string;
using std::array;
using std::vector;
/*
namespace std {
    template<>
    class hash<tuple<int,int>> {
    public:
        std::size_t operator()(const tuple<int,int>& t) const {
            return get<0>(t) * 10 + get<1>(t);
        }
    };
}
*/

template< typename T1, typename T2=std::nullptr_t >
struct Std{
typedef unordered_set<T1, std::hash<T1>> Set;
typedef unordered_map<T1, T2> Map;
};

template< typename T1, typename T2=std::nullptr_t >
struct Google{
    typedef dense_hash_set<T1, std::hash<T1>> Set;
    typedef dense_hash_map<T1, T2, std::hash<T1> > Map;
};




int main(int argc, char** argv) {
/*
typedef unordered_map<int, unordered_set<string>> NodeMapT;
NodeMapT b;
b[3] = unordered_set<string>();
b[3].insert("Hello");

 auto t1 = make_tuple("hi", "you");
 auto t2 = make_tuple("hi", "man");

 std::cout << "t1 =?= t2 " << (t1 == t2) << "\n";

 //unordered_map<tuple<int,int>, unordered_map<tuple<int,int>,string> > tm;

 //tm[ make_tuple(0,0) ] = unordered_map<tuple<int,int>,string>();
 //tm[ make_tuple(0,0) ][ make_tuple(1,1) ] = "b+";

 // unordered_map<tuple<int,int>, string> ot = { { make_tuple(0,0) , "b+" } };

 vector<int> vec = {1,2,3,4};

 //  std::cout << MultiOpt::flipDict[ make_tuple(0,0) ][ make_tuple(0,0) ]<< "\n";
 // std::cout << tm[ make_tuple(0,0) ][ make_tuple(1,1) ] << "\n";
 // std::cout << ot[ make_tuple(0,0) ] << "\n";
 // std::cout << v[0] << "\n";
 using std::vector;
 using std::unordered_set;
 using boost::bimap;

 using boost::bimaps::vector_of_relation;
 using boost::bimaps::vector_of;
 using boost::bimaps::unordered_set_of;
 using std::hash;

 using std::make_pair;

 typedef bimap< vector_of< size_t >, vector_of<string>, vector_of_relation > VertexBiMapT;
 //        vector_of< size_t >,
 //       vector_of< string >
 //       > VertexBiMapT;

    typedef VertexBiMapT::value_type VertexBiKey;

    VertexBiMapT vBiMap;

    vBiMap.push_back( VertexBiKey(0, "hi") );
    vBiMap.push_back( VertexBiKey(1, "hi1") );
    vBiMap.push_back( VertexBiKey(2, "hi2") );
    for( auto e : vBiMap.left) {
        std::cout << "left : "<< e.first << ", right : " << e.second << "\n";
    }


    //VertexBiMapT::right_value_type t = VertexBiMapT::right_value_type("WTF");
    //std::cout << t << "\n";

    using boost::adjacency_list;
    using boost::vecS;
    using boost::setS;
    using boost::directedS;
    using boost::undirectedS;
    using boost::add_edge;
    using boost::labeled_graph;

class Vertex{
public:
string name;
size_t idx;
};

typedef labeled_graph<
    adjacency_list<vecS, vecS, undirectedS, Vertex>,
    size_t
    > graphT;
graphT G;

auto v = add_vertex(1,G); G[v].name = "one"; G[v].idx = 1;
v = add_vertex(8,G);  G[v].name = "eight"; G[v].idx = 8;
  v = add_vertex(8,G);  G[v].name = "eight"; G[v].idx = 8;
v = add_vertex(8,G);  G[v].name = "eight"; G[v].idx = 8;
v = add_vertex(15,G);  G[v].name = "fifteen"; G[v].idx = 15;

add_edge_by_label(1,8,G);
add_edge_by_label(8,8,G);
add_edge_by_label(8,15,G);

auto ve = boost::vertices(G);
std::cout << "VERTICES : \n";
    for ( auto it = ve.first; it != ve.second; ++it ) {
        std::cout << *it << " ";
    } std::cout << "\n";

auto av = boost::adjacent_vertices(8,G.graph());
std::cout << "ADJ:" << "\n";
    for( auto it = av.first; it != av.second; ++it ) {
        std::cout << *it << " ";
    }
*/

    /*
cerr << "using std::set " << "\n";
seed_val = 1234;
rng.seed(seed_val);
std::uniform_int_distribution<unsigned int> uint_dist;
std::vector< std::set< unsigned int > > sets;

auto time1 = std::chrono::high_resolution_clock::now();

cerr << "creating sets\n";
for( size_t i = 0; i < 200; ++i ) {

sets.push_back( std::set<unsigned int>() );

for( size_t j = 0; j < 1000; ++j ) {
sets.back().insert( uint_dist(rng) );
}

}
cerr << "done\n";

cerr << "comparing sets\n";

for ( auto iit = sets.begin(); iit != sets.end(); ++iit ) {
auto jit = iit;
jit++;
for ( ; jit != sets.end(); ++jit) {
vector<unsigned int> res;
std::set_union( iit->begin(), iit->end(), jit->begin(), jit->end(), std::back_inserter(res));
//cerr << "size of union is " << res.size() << " ";
}
}
auto time2 = std::chrono::high_resolution_clock::now();
std::cout << "f() took "
<< std::chrono::duration_cast<std::chrono::milliseconds>(time2-time1).count() / 1000.0
<< " milliseconds\n";

cerr << "done\n";





cerr << "using std::unordered_set " << "\n";

time1 = std::chrono::high_resolution_clock::now();

 std::vector< std::unordered_set< unsigned int > > hsets;

cerr << "creating sets\n";
for( size_t i = 0; i < 200; ++i ) {

hsets.push_back( std::unordered_set<unsigned int>() );

for( size_t j = 0; j < 1000; ++j ) {
hsets.back().insert( uint_dist(rng) );
}

}
cerr << "done\n";

cerr << "comparing sets\n";

for ( auto iit = hsets.begin(); iit != hsets.end(); ++iit ) {
auto jit = iit;
jit++;
for ( ; jit != hsets.end(); ++jit) {
  unordered_set<unsigned int> res;
  for ( auto e : *iit ) { res.insert(e); }
  for ( auto e : *jit ) { res.insert(e); }
//std::set_union( iit->begin(), iit->end(), jit->begin(), jit->end(), std::back_inserter(res));
//cerr << "size of union is " << res.size() << " ";
}
}

time2 = std::chrono::high_resolution_clock::now();
std::cout << "f() took "
<< std::chrono::duration_cast<std::chrono::milliseconds>(time2-time1).count() / 1000.0
<< " milliseconds\n";
*/



    std::set<int> a = {1,2,3,4,5};
    std::set<int> b = {3,4,7,8,9};

    std::set<int> c = {1,7,12};

    std::set_union( a.begin(), a.end(), a.begin(), a.end(), std::inserter(c, c.begin()) );

    cerr << "union : ";
    for ( auto e : c ) {
        cerr << e <<  " ";
    } cerr << "\n";


    std::vector<int> d = {1};//,2,3,4,5,6,7};

    for( auto it = d.begin(); it < d.end(); it+=2 ) {
        if ( std::distance(it, d.end()) > 1 ) {
            cerr << *it << ", " << *(it+1) << " ";
        } else {
            cerr << *it << " ";
        }
    } cerr << "\n";




    Std<int>::Set h;
    Google<int, std::vector<Derivation> >::Map m;
    m.set_empty_key(-1);

    for ( size_t i = 0; i < 20; ++i) {
        h.insert(i);
        m[i] = std::vector<Derivation>();
        m[i].push_back( Derivation() );
    }

    for( auto& e : h ) {
        cerr << e << " " ;
    } cerr << "\n";

    for ( auto& e : m ) {
        cerr << e.first << ", " << e.second.front() << "\n";
    }
    cerr << "\n";
    return 0;
}
