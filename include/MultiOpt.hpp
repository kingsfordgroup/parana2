#ifndef MULTI_OPT_HPP
#define MULTI_OPT_HPP

#include <fstream>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include <Bpp/Phyl/Tree.h>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/functional/hash.hpp>
#include <cln/rational.h>
#include <cln/integer.h>
#include <cln/float.h>
#include <cln/real.h>
#include <cln/rational_io.h>
#include <cln/integer_io.h>
#include <cln/float_io.h>
#include <limits>
#include <sstream>
#include <exception>
#include <sstream>
#include "mpreal.h"
#include "TreeUtils.hpp"
#include "CountedDerivation.hpp"
#include "CostClass.hpp"
#include "Derivation.hpp"
#include "FHGraph.hpp"
#include "GraphUtils.hpp"

/** Google's dense hash set and hash map **/
#include <google/dense_hash_set>
#include <google/dense_hash_map>
#include <google/sparse_hash_map>

namespace std {
    template<>
    class hash<tuple<int,int>> {
    public:
        std::size_t operator()(const tuple<int,int>& t) const {
            size_t seed = 0;
            boost::hash_combine(seed, get<0>(t));
            boost::hash_combine(seed, get<1>(t));
            return seed;
        }
    };
}

namespace std {
    template<>
    class hash<tuple<size_t,size_t>> {
    public:
        std::size_t operator()(const tuple<int,int>& t) const {
            size_t seed = 0;
            boost::hash_combine(seed, get<0>(t));
            boost::hash_combine(seed, get<1>(t));
            return seed;
        }
    };
}


namespace std {
    template<>
    class hash<tuple<bool,bool>> {
    public:
        std::size_t operator()(const tuple<bool,bool>& t) const {
            size_t seed = 0;
            boost::hash_combine(seed, get<0>(t));
            boost::hash_combine(seed, get<1>(t));
            return seed;
        }
    };
}

namespace MultiOpt {

    using cln::cl_I;
    using cln::cl_F;
    using cln::cl_RA;
    using cln::cl_R;
    using cln::cl_float;
    using google::dense_hash_set;
    using google::dense_hash_map;
    using google::sparse_hash_map;
    using std::map;
    using std::cout;
    using std::unordered_map;
    using std::unordered_set;
    using std::tuple;
    using std::tie;
    using std::make_tuple;
    using std::string;
    using std::get;
    using bpp::Tree;
    using Utils::TreeInfo;
    using Utils::Trees::differentExtantNetworks;
    using boost::assign::map_list_of;
    using std::unique_ptr;
    using std::ostringstream;
    using mpfr::mpreal;
    using Utils::Trees::TreePtrT;

    template< typename T1, typename T2=std::nullptr_t >
    struct Std{
        typedef unordered_set<T1, std::hash<T1>> Set;
        typedef unordered_map<T1, T2, std::hash<T1>> Map;
    };

    template< typename T1, typename T2=std::nullptr_t >
    struct Google{
        typedef dense_hash_set<T1, std::hash<T1>> Set;
        typedef dense_hash_map<T1, T2, std::hash<T1>> Map;
    };

    typedef tuple<double, vector<size_t> > dvsT;
    typedef tuple<double, size_t, vector<size_t> > edvsT;

    template <typename elemT>
    class QueueCmp {
    public:
        bool operator() ( const elemT& lhs, const elemT& rhs ) const  {
            return get<0>(lhs) > get<0>(rhs);
        }
        bool operator() ( const elemT& lhs, const elemT& rhs ) {
            return get<0>(lhs) > get<0>(rhs);
        }
    };

    typedef tuple<double, vector<size_t> > dvsT;
    typedef unordered_map< size_t, unordered_map<size_t, Derivation>> slnDictT;
    typedef Google< size_t, vector< tuple<double, cl_I> > >::Map countDictT;

    typedef tuple<bool,bool> flipTupleT;
    typedef tuple<double, string> costRepT;
    typedef unordered_map< flipTupleT, unordered_map<flipTupleT, string> > flipMapT;
    typedef unordered_map< flipTupleT, unordered_map<flipTupleT, costRepT> > costMapT;
    typedef unordered_map< flipTupleT, unordered_map<flipTupleT, std::function< costRepT (const double&, const double&) > > > costMapFunT;
    typedef unordered_map< flipTupleT, unordered_map<bool, costRepT> > selfCostMapT;
    typedef unordered_map< flipTupleT, unordered_map<bool, std::function<costRepT (const double&) > > > selfCostMapFunT;

    //    Google< tuple<size_t,size_t> >::Set recursedStore;
    //    Google<size_t, vector<Derivation>*>::Map cand;


    costMapT getCostDict ( double cc, double dc, bool directed );

    costMapFunT getCostFunDict ( double cc, double dc, bool directed );

    selfCostMapT getSelfLoopCostDict ( double cc, double dc, bool directed );

    selfCostMapFunT getSelfLoopCostFunDict ( double cc, double dc, bool directed );

    string flipType( const FlipKey& hvert, const FlipKey& tvert );

    template<typename GT>
    unordered_set<int> projectToReversedGraph( unique_ptr<ForwardHypergraph>& H, GT& G );

    void topologicalOrder( unique_ptr<ForwardHypergraph>& H, vector<size_t>& order );

    /**
     *  Compute the penalty for this edge to exist based on difference
     *  between the existence intervals of the endpoints and the
     *  penalty factor.
     */
    template <typename T>
    double existencePenalty( const TreeInfo& ti, const T& vert, const double& penalty, const double& travWeight );

    unique_ptr<ForwardHypergraph>  buildSolutionSpaceGraph( const TreePtrT& t,
                                                            const TreeInfo& ti,
                                                            double cc,
                                                            double dc,
                                                            double penalty,
                                                            bool directed );


    template< typename GT >
    void leafCostDict( unique_ptr<ForwardHypergraph>& H, TreePtrT& T, GT& G, bool directed, double cc, double dc, slnDictT& slnDict );

    tuple<double, cl_I> getCostCount( const unordered_map<size_t, vector<CostClass> >& tkd,
                                      const vector<size_t>& bp,
                                      const size_t& eid,
                                      unique_ptr<ForwardHypergraph>& H );

    double getCost( const unordered_map<size_t, vector<CostClass> >& tkd,
                    const vector<size_t>& bp,
                    const size_t& eid,
                    unique_ptr<ForwardHypergraph>& H );

    cl_I getCount( const unordered_map<size_t, vector<CostClass> >& tkd,
                   const vector<size_t>& bp,
                   const size_t& eid,
                   unique_ptr<ForwardHypergraph>& H );

    vector<CostClass> computeKBest(const size_t& vid,
                                   const size_t& k,
                                   const unordered_map<size_t, vector<CostClass> >& tkd,
                                   unique_ptr<ForwardHypergraph>& H);


    // 0 : 20 23 25
    // 1 : 15 16 34
    // 2 : 10 12 13
    // rank, #
    // (0, 3000)
    //
    vector< tuple<double, cl_I> > countEdgeSolutions( const double& ecost,
                                                      const vector<size_t>& tailNodes,
                                                      countDictT& countDict,
                                                      const size_t& k,
                                                      bool printMe,
                                                      unique_ptr<ForwardHypergraph>& H,
                                                      TreePtrT& t
                                                       );

    double estimateThermodynamicBeta( const vector<tuple<double, cl_I>>& slnVecIn, const double& emin );

    vector<double> getAlphas( const vector< tuple<double, cl_I> >& slnVec, const cl_I& numPaths );

    void insideOutside( unique_ptr<ForwardHypergraph>& H, TreePtrT& t, TreeInfo& ti, double penalty, const vector<size_t>& order, slnDictT& slnDict, countDictT& countDict, const string& outputName );

    FlipKey keyForAction( const FlipKey& fk , const string& ft );

    void viterbiCount( unique_ptr<ForwardHypergraph>& H,
                       TreePtrT& t,
                       TreeInfo& ti,
                       double penalty,
                       const vector<size_t>& order,
                       slnDictT& slnDict,
                       countDictT& countDict,
                       const size_t& k,
                       const string& outputName,
                       const vector<FlipKey>& outputKeys );

    void viterbiCountNew( unique_ptr<ForwardHypergraph>& H, TreePtrT& t, TreeInfo& ti, double penalty, const vector<size_t>& order,
                          slnDictT& slnDict, countDictT& countDict, const size_t& k,
                          const string& outputName, const vector<FlipKey>& outputKeys );

    void viterbi( unique_ptr<ForwardHypergraph>& H, TreePtrT& t, TreeInfo& ti, double penalty, const vector<size_t>& order, slnDictT& slnDict );

    /** Alg 3 from the paper */
    typedef std::pair<size_t, Derivation> TaggedDerivT;

    //Google< tuple<size_t,size_t> >::Set recursedStore;
    //Google<size_t, vector<Derivation>*>::Map cand;
    //unordered_map< size_t, vector<Derivation>* > cand;

    void initKBest();

    template<typename T>
    void printVector( const vector<T>& v );

    void getCandidates( const unique_ptr<ForwardHypergraph>& H, TreeInfo& ti, double penalty, slnDictT& d, size_t v, size_t k, vector<Derivation>* kbest );

    bool sameEffect( const Derivation& d0, const Derivation& d1 );

    void lazyKthBest( unique_ptr<ForwardHypergraph>& H, TreeInfo& ti, double penalty, slnDictT& d, size_t v, size_t k, size_t kp );

    void lazyNext( unique_ptr<ForwardHypergraph>& H, TreeInfo& ti, double penalty, slnDictT& d, vector<Derivation>* locCand, size_t eind,
                   const vector<size_t>& j, size_t kp);

    void lazyKthBest( unique_ptr<ForwardHypergraph>& H, TreeInfo& ti, double penalty, slnDictT& d, size_t v, size_t k, size_t kp );
}




#endif // MULTI_OPT_HPP
