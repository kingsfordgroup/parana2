#ifndef TREEUTILS_HPP
#define TREEUTILS_HPP

/** Standard Includes */
#include <map>
#include <set>
#include <tuple>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <memory>
#include <tuple>
#include <vector>
//#include "FlipKey.hpp"

/** Bio++ **/
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Io/NexusIOTree.h>
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Exceptions.h>

namespace Utils {
    using std::get;
    using std::cerr;
    using std::cout;
    using std::string;
    using std::tuple;
    using std::make_tuple;
    using std::unordered_map;
    using std::unordered_set;
    using std::endl;
    using std::vector;
    using std::unique_ptr;
    using bpp::Newick;
    using bpp::Tree;
    using bpp::Exception;

    class TreeInfo {
    public:
        typedef unordered_map<int, unordered_set<string>> NodeMapT;
        typedef unordered_map<int, unordered_set<int>> NodeIndMapT;
        typedef unordered_map<int, tuple<double,double>> NodeIntervalMapT;
        NodeIndMapT subnodes;
        NodeIndMapT leaves;
        NodeMapT enets;
        NodeIntervalMapT extantInterval;
        string name;

        TreeInfo( const string n ) : name(n), subnodes(NodeIndMapT()),
                                     leaves(NodeIndMapT()),
                                     enets(NodeMapT()),
                                     extantInterval(NodeIntervalMapT()) { cout << "creating tree info: " << name << "\n"; }

        bool inSubnodesOf( int u,  int v ) {
            return subnodes[u].find(v) != subnodes[u].end();
        }
        bool inSubnodesOf( int u,  int v ) const {
            return subnodes.find(u) != subnodes.cend() && subnodes.find(u)->second.find(v) != subnodes.find(u)->second.cend();
        }

        double intervalDistance( int u, int v ) {
            //return static_cast<const TreeInfo&>(*this).intervalDistance(u,v);
            auto birthU = get<0>( extantInterval[u] );
            auto deathU = get<1>( extantInterval[u] );
            auto birthV = get<0>( extantInterval[v] );
            auto deathV = get<1>( extantInterval[v] );
            // If the ranges overlap, the distance is 0
            if ( (birthU <= deathV) && (birthV <= deathU) ) { return 0.0; }
            // otherwise
            if ( birthU > deathV ) {
                return birthU - deathV;
            }

            if ( birthV > deathU ) {
                return birthV - deathU;
            }
        }

        double intervalDistance( int u, int v ) const {
            return const_cast<TreeInfo&>(*this).intervalDistance(u,v);
        }
    };

    bool advanceElems( vector<size_t>& ptrs, const vector<size_t>& sizes );

    template <typename pqT, typename pqCompT>
    bool appendNext( double score,
                     const vector<size_t>& inds,
                     const vector<size_t>& sizes,
                     vector<pqT>& pq,
                     pqCompT& pqComp,
                     std::function< double(const vector<size_t>&) >& computeScore );

    namespace Trees {
        typedef unique_ptr<bpp::TreeTemplate<bpp::Node>> TreePtrT;

        bool differentExtantNetworks( const TreeInfo& tinfo, int u, int v );

        void labelTree( const TreePtrT& tree );

        TreePtrT readNewickTree( const std::string& treeName );

        const std::string getName( TreePtrT& t, int nid);

        string getExtantNetwork(const string& s);

        void prepareTree( TreePtrT& t, TreeInfo& ti, int nid );
    }

}

#endif // TREEUTILS_HPP
