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

/** Bio++ **/
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Io/NexusIoTree.h>
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Exceptions.h>

//#include "FlipKey.hpp"
#include "ParanaCommon.hpp"

// forward declaration
class FlipKey;

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
    using std::shared_ptr;
    using bpp::Newick;
    using bpp::Tree;
    using bpp::Exception;

    double round3( double num );

    class ExistenceInterval {
        public:
            double birth;
            double death;
            ExistenceInterval( ) : birth(std::numeric_limits<double>::infinity()), 
                                   death(-std::numeric_limits<double>::infinity()) {}

            ExistenceInterval( double _birth, double _death ) :
                               birth(_birth), death(_death) {}

            double distance( const ExistenceInterval& other ) {
                if ( birth <= other.death and other.birth <= death ) {
                    return 0.0;
                }
                if ( birth > other.death ) {
                    return birth - other.death;
                }
                if ( other.birth > death ) {
                    return other.birth - death;
                }
            }
    };

    class TreeInfo {
    public:
        typedef shared_ptr<bpp::TreeTemplate<bpp::Node>> TreePtrT;
        typedef unordered_map<int, unordered_set<string>> NodeMapT;
        typedef unordered_map<int, unordered_set<int>> NodeIndMapT;
        typedef unordered_map<int, ExistenceInterval> NodeIntervalMapT;
        TreePtrT tree;
        NodeIndMapT subnodes;
        NodeIndMapT leaves;
        NodeMapT enets;
        NodeMapT subspec;
        NodeIntervalMapT extantInterval;
        string name;



       TreeInfo( const string n, TreePtrT t ) : name(n), tree(t), subnodes(NodeIndMapT()),
                                     leaves(NodeIndMapT()),
                                     enets(NodeMapT()),
                                     extantInterval(NodeIntervalMapT()) { }

        bool inSubnodesOf( int u,  int v ) {
            return subnodes[u].find(v) != subnodes[u].end();
        }
        bool inSubnodesOf( int u,  int v ) const {
            return subnodes.find(u) != subnodes.cend() && subnodes.find(u)->second.find(v) != subnodes.find(u)->second.cend();
        }

        double intervalDistance( int u, int v ) {
            return extantInterval[u].distance(extantInterval[v]);
        }

        double intervalDistance( int u, int v ) const {
            return const_cast<TreeInfo&>(*this).intervalDistance(u,v);
        }
    };

    bool advanceElems( vector<size_t>& ptrs, const vector<size_t>& sizes );

    template <typename pqT, typename pqCompT>
    bool appendNext( double score,
                     const vector<size_t, StackAllocator<size_t>>& inds,
                     const vector<size_t, StackAllocator<size_t>>& sizes,
                     vector<pqT>& pq,
                     pqCompT& pqComp,
                     std::function< double(const vector<size_t, StackAllocator<size_t>>&) >& computeScore );

    template <typename pqT>
    bool appendNextWithEdge( const size_t& eid,
                             const vector<size_t, StackAllocator<size_t>>& inds,
                             const vector<size_t, StackAllocator<size_t>>& sizes,
                             pqT& pq,
                             std::function< double(const size_t& eid, const vector<size_t, StackAllocator<size_t>>&) >& computeScore );
    /*
    template <typename pqT>
    bool appendNextWithEdgeOrig( const size_t& eid,
                                 const vector<size_t>& inds,
                                 const vector<size_t>& sizes,
                                 pqT& pq,
                                 std::function< double(const size_t& eid, const vector<size_t>&) >& computeScore );
    */

    namespace Trees {
        typedef shared_ptr<bpp::TreeTemplate<bpp::Node>> TreePtrT;

        bool differentExtantNetworks( const TreeInfo& tinfo, int u, int v );

        void labelTree( TreePtrT& tree );

        bool sameSpecies(const TreePtrT& t, int u, int v);

        bool isSpeciationEvent(const FlipKey& k, const TreeInfo& ti);

        std::tuple<int, int, int, int> correspondingSpecies(const FlipKey& k, const TreeInfo& ti);

        bool isDescendantSpecies(const TreeInfo& ti, int u, int v);

        TreePtrT readNewickTree( const std::string& treeName );

        const std::string getName( const TreePtrT& t, int nid);

        string getExtantNetwork(const string& s);

        void prepareTree( TreePtrT& t, TreeInfo& ti, int nid );
    }

    string stringForKey (const FlipKey& key, const Trees::TreePtrT& t);
}

#endif // TREEUTILS_HPP
