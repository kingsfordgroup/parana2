#include "MultiOpt.hpp"
#include "TreeUtils.hpp"
#include <boost/heap/fibonacci_heap.hpp>
#include <boost/heap/pairing_heap.hpp>
#include <boost/heap/binomial_heap.hpp>
#include <boost/heap/skew_heap.hpp>
#include <Bpp/Phyl/Io/Nhx.h>
#include "ParanaCommon.hpp"
#include "FlipKey.hpp"

using boost::heap::fibonacci_heap;
using boost::heap::pairing_heap;

namespace Utils {

    double round3( double num ) {
        double result = num * 100000;
        result = std::floor(result);
        result = result / 100000;
        return result;
    }

    
    string stringForKey (const FlipKey& key, const Trees::TreePtrT& t) {
        auto uname = t->getNodeName(key.u());
        auto vname = t->getNodeName(key.v());
        if (uname > vname) {
            auto tmp = vname;
            vname = uname;
            uname = tmp;
        }

        string fstr = "";
        if ( key.state() == FlipState::both ) {
            fstr = "<-->";
        } else if ( key.state() == FlipState::none ) {
            fstr = "X";
        } else {
            std::abort();
        }
        return "[" + uname + ", " + vname + "] : " + fstr;
    } 
    
   
    bool advanceElems( vector<size_t>& ptrs, const vector<size_t>& sizes ) {

        assert( ptrs.size() == sizes.size() );
        // Find the last element
        size_t ctr = sizes.size()-1;
        for ( auto rp = ptrs.rbegin(); rp != ptrs.rend(); ++rp ) {
            if ( *rp < sizes[ctr]-1 ) {
                ++(*rp);
                return true;
            } else {
                if ( rp+1 == ptrs.rend() ) {
                    return false;
                } else {
                    (*rp) = 0;
                }
            }
            --ctr;
        }

    }

    template <typename pqT, typename pqCompT>
    bool appendNext( double score,
                     const vector<size_t, StackAllocator<size_t>>& inds,
                     const vector<size_t, StackAllocator<size_t>>& sizes,
                     vector<pqT>& pq,
                     pqCompT& pqComp,
                     std::function< double(const vector<size_t, StackAllocator<size_t>>&) >& computeScore ) {

        size_t i = 0;
        while ( i < inds.size() ) {
            vector<size_t, StackAllocator<size_t>> newInds(inds);
            //newInds.reserve(inds.size());
            //std::copy(inds.begin(), inds.end(), newInds.begin());
            newInds[i]++;
            if ( newInds[i] < sizes[i] ) {
                pq.push_back( make_tuple( computeScore(newInds), newInds ) );
                std::push_heap( pq.begin(), pq.end(), pqComp );
            }
            if (inds[i] != 0) { return true; }
            i += 1;
        }

        return true;
    }

    template <typename pqT>//, typename pqCompT>
    bool appendNextWithEdge( const size_t& eid,
                             const vector<size_t, StackAllocator<size_t>>& inds,
                             const vector<size_t, StackAllocator<size_t>>& sizes,
                             pqT& pq,
                             std::function< double(const size_t& eid, const vector<size_t,StackAllocator<size_t>>&) >& computeScore ) {
        size_t i = 0;
        size_t ninds = inds.size();
        while ( i < ninds ) {
            vector<size_t, StackAllocator<size_t>> newInds(inds);
            newInds[i]++;
            if ( newInds[i] < sizes[i] ) {
                // Custom fix -- check to see if Boost folks fix emplace operation
                // pq.emplace( computeScore(eid, newInds), eid, newInds );
                pq.emplace( MultiOpt::EdgeDerivation(computeScore(eid, newInds), eid, newInds) );
            }
            if (inds[i] != 0) { return true; }
            i += 1;
        }

        return true;
    }

    /*
    template <typename pqT>//, typename pqCompT>
    bool appendNextWithEdgeOrig( const size_t& eid,
                                 const vector<size_t>& inds,
                                 const vector<size_t>& sizes,
                                 pqT& pq,
                                 std::function< double(const size_t& eid, const vector<size_t>&) >& computeScore ) {

        size_t i = 0;
        while ( i < inds.size() ) {
            vector<size_t> newInds(inds);
            newInds[i]++;
            if ( newInds[i] < sizes[i] ) {
                double score = computeScore(eid, newInds);
                bool exists = false;
                auto it = pq.ordered_begin();
                while( it != pq.ordered_end() && !exists ) {
                    double dscore = get<0>(*it);
                    size_t oeid = get<1>(*it);
                    if ( dscore == score && eid == oeid ) {
                        auto oinds = get<2>(*it);
                        size_t sct = 0;
                        for( size_t j = 0; j < oinds.size(); ++j ) {
                            if( newInds[j] != oinds[j] ){ break; } else { sct += 1; }
                        }
                        if ( sct == oinds.size() ) { exists = true; }
                    }
                    if ( dscore > score ) { break; }
                    ++it;
                }

                if (!exists) {
                    pq.push( make_tuple( computeScore(eid, newInds), eid, newInds ) );
                }
            }
            //if (inds[i] != 0) { return true; }
            i += 1;
        }

        return true;
    }
*/

    typedef tuple<double, vector<size_t> > dvsT;

    class QueueCmp {
    public:
        bool cmpArrays( const vector<size_t>& a, const vector<size_t>& b) {
            assert(a.size() == b.size());
            size_t i = 0;
            while (i < a.size() && a[i] == b[i] ) {
                i += 1;
            }
            if (i == a.size()) {
                return true;
            } else {
                return a[i] > b[i];
            }
        }

        bool operator() ( const dvsT& lhs, const dvsT& rhs ) {
                return get<0>(lhs) >= get<0>(rhs);
        }
    };



    namespace Trees {
        typedef shared_ptr<bpp::TreeTemplate<bpp::Node>> TreePtrT;


        bool differentExtantNetworks( const TreeInfo& tinfo, int u, int v ) {
            //auto specU = dynamic_cast< bpp::BppString* >(tinfo.tree->getNodeProperty(u, "S"))->toSTL();
            //auto specV = dynamic_cast< bpp::BppString* >(tinfo.tree->getNodeProperty(v, "S"))->toSTL();
            //return specU != specV;
            if ( u == v ) { return false; }
            auto& enetsU = tinfo.enets.find(u)->second;
            auto& enetsV = tinfo.enets.find(v)->second;
            
            for (auto uit = enetsU.cbegin(); uit != enetsU.cend(); ++uit) {
                if ( enetsV.find(*uit) != enetsV.cend() ) { return false; }
            }
            return true;
        }

        /**
         * Returns `true` if the event occuring at node `k` is a speciation event
         * (i.e. both k.u() and k.v() have decendants in different species) and `false`
         * otherwise.
         * @param  k The hypervertex to check
         * @param  t The phylogenetic tree 
         * @return   `true` if the decendants of this hypervertex represent
         *           a speciation and `false` otherwise
         */
        bool isSpeciationEvent(const FlipKey& k, const TreeInfo& ti) {
            auto u = k.u();
            auto v = k.v();
            const TreePtrT& t = ti.tree;
            // If either node is a leaf, they can't have decendant species
            if (t->isLeaf(u) or t->isLeaf(v)) { return false; }

            // Get the child nodes
            int lu, ru, lv, rv;
            lu = ru = lv = rv = -1;            
    
            auto sons = t->getSonsId(u);
            lu = sons[0]; ru = sons[1]; 
            if (lu > ru) { std::swap(lu, ru); }

            sons = t->getSonsId(v);
            lv = sons[0]; rv = sons[1]; 
            if (lv > rv) { std::swap(lv, rv); }

            auto specu = dynamic_cast< bpp::BppString* >(ti.tree->getNodeProperty(u, "S"))->toSTL();
            auto specv = dynamic_cast< bpp::BppString* >(ti.tree->getNodeProperty(v, "S"))->toSTL();
            auto speclu = dynamic_cast< bpp::BppString* >(ti.tree->getNodeProperty(lu, "S"))->toSTL();
            auto speclv = dynamic_cast< bpp::BppString* >(ti.tree->getNodeProperty(lv, "S"))->toSTL();
            auto specru = dynamic_cast< bpp::BppString* >(ti.tree->getNodeProperty(ru, "S"))->toSTL();
            auto specrv = dynamic_cast< bpp::BppString* >(ti.tree->getNodeProperty(rv, "S"))->toSTL();                                    

            bool parentsAreSameSpecies = specu == specv;
            bool uSpeciates = speclu != specru;
            bool vSpeciates = speclv != specrv;
            bool speciesMatch = ((speclu == speclv) or (speclu == specrv)) and
                                ((specru == speclv) or (specru == specrv)); 


            bool speciationEvent = (parentsAreSameSpecies and uSpeciates and vSpeciates and speciesMatch);

            return speciationEvent;
        }

        /**
         * NOTE: Assumes isSpeciationEvent(k, ti) !!
         * For the node k = ({u, v}, s), returns the children of u and v corresponding
         * to the first species as the first two elements of the tuple and the children
         * corresponding to the second species as the last two elements of the tuple.
         * @param  k  Hypervertex whose children are to be put into correspondence
         * @param  ti The tree info
         */
        std::tuple<int, int, int, int> correspondingSpecies(const FlipKey& k, const TreeInfo& ti) {
            auto u = k.u();
            auto v = k.v();
            const TreePtrT& t = ti.tree;

            // Get the child nodes
            int lu, ru, lv, rv;
            lu = ru = lv = rv = -1;            
    
            auto sons = t->getSonsId(u);
            lu = sons[0]; ru = sons[1]; 
            if (lu > ru) { std::swap(lu, ru); }

            sons = t->getSonsId(v);
            lv = sons[0]; rv = sons[1]; 
            if (lv > rv) { std::swap(lv, rv); }

            auto speclu = dynamic_cast< bpp::BppString* >(ti.tree->getNodeProperty(lu, "S"))->toSTL();
            auto speclv = dynamic_cast< bpp::BppString* >(ti.tree->getNodeProperty(lv, "S"))->toSTL();
            auto specru = dynamic_cast< bpp::BppString* >(ti.tree->getNodeProperty(ru, "S"))->toSTL();
            auto specrv = dynamic_cast< bpp::BppString* >(ti.tree->getNodeProperty(rv, "S"))->toSTL();                                    

            if (speclu == speclv) {
                //std::cerr << "returning (" << speclu << ", " << speclv <<
                //             ", " << specru << ", " << specrv << ")\n";
                return make_tuple(lu, lv, ru, rv);
            } else {
                //std::cerr << "returning (" << speclu << ", " << specrv <<
                //             ", " << specru << ", " << speclv << ")\n";
                return make_tuple(lu, rv, ru, lv);
            }
        }

        TreePtrT readNewickTree( const std::string& treeName ) {
            unique_ptr<bpp::Nhx> newickReader( new bpp::Nhx(true) ); //No comment allowed!
            //newickReader->enableExtendedBootstrapProperty("name");
            try {
                TreePtrT tree( newickReader->read(treeName) ); // Tree in file
                for ( auto pn : tree->getNodePropertyNames( tree->getRootId()) ) {
                    std::cerr << "property : " << pn << "\t";
                    std::cerr << "value : " << dynamic_cast<bpp::BppString*>(tree->getNodeProperty(tree->getRootId(), pn))->toSTL() << "\n";
                }

                return tree;
            } catch(std::exception& e) {
                cerr << "Caught Exception: [" << e.what() << "]\n";
                abort();
            }
        }


        const std::string getName( const TreePtrT& t, int nid) {
            if ( t->hasNodeProperty(nid,"GN") ) {
                auto name = dynamic_cast<bpp::BppString*>( t->getNodeProperty(nid,"GN"))->toSTL();
                return name;
            } else if (t->hasNodeName(nid)) {
                return t->getNodeName(nid);
            }
        }

        void labelTree( TreePtrT& tree ) {
            for ( auto nid : tree->getNodesId() ) {
                tree->setNodeName(nid, getName(tree, nid) );
            }
        }

        string getExtantNetwork(const string& s) {
            return ( s.find("LOST") != string::npos ) ? "LOST" : s.substr( s.find_last_of('_') );
        }

        string getSpeciesName(const TreePtrT& t, int nid) {
            return dynamic_cast<bpp::BppString*>( t->getNodeProperty(nid,"S") )->toSTL();
        } 

        bool isDescendantSpecies(const TreeInfo& ti, int u, int v)  {

            if ( ti.tree->isLeaf(v) ) { return false; }
            auto sonIDs = ti.tree->getSonsId(v);

            auto& subspecLV = ti.subspec.find(sonIDs[0])->second;
            auto& subspecRV = ti.subspec.find(sonIDs[1])->second;
            
            auto specU = getSpeciesName(ti.tree, u);

            auto res = ( (subspecLV.find( specU ) != subspecLV.end()) or 
                         (subspecRV.find( specU ) != subspecRV.end()) );
            return res;

            /*
            auto& enetsU = ti.enets.find(u)->second;
            auto& enetsV = ti.enets.find(v)->second;

            // If u has more descendant species than v, then it can't be a descendant of v
            if (enetsU.size() > enetsV.size()) { return false; }

            // Every species below u must be represented in v
            for (auto uit = enetsU.cbegin(); uit != enetsU.cend(); ++uit) {
                if ( enetsV.find(*uit) == enetsV.cend() ) { return false; }
            }
            return true;
            */
        }

        bool sameSpecies(const TreePtrT& t, int u, int v) {
            return  getSpeciesName(t, u) == getSpeciesName(t, v);
        }

        void prepareTree( TreePtrT& t, TreeInfo& ti, int nid ) {
            // if the current node is not the root
            if ( nid != t->getRootId() ) {
                
                auto fid = t->getFatherId(nid);
                
                if ( ti.extantInterval.find(fid) == ti.extantInterval.end() ) {
                    std::cerr << "FATAL: Must know parent interval before computing "
                              << "child interval\n";
                    std::abort();
                }

                auto parentDeathT = ti.extantInterval[fid].death;
                double fdist = 0.0;
                if ( t->isLeaf(nid) ) {
                    fdist = std::numeric_limits<double>::infinity();
                } else if ( t->hasDistanceToFather(nid) ) {
                    fdist = t->getDistanceToFather(nid);
                }
                ti.extantInterval[ nid ] = ExistenceInterval( parentDeathT, parentDeathT + fdist );
            }

            if (t->isLeaf(nid)) {
                ti.leaves[nid] = {nid};
                ti.subnodes[nid] = {nid};
                string enet = getSpeciesName(t,nid); //getExtantNetwork(getName(t,nid));

                // Skip lost nodes
                ti.enets[nid] = {enet};
                ti.subspec[nid] = {enet};
                
                // Consider lost nodes
                //ti.enets[nid] = {enet};
            } else {
                ti.leaves[nid] = {};
                ti.subnodes[nid] = { nid };
                ti.enets[nid] = unordered_set<string>();
                ti.subspec[nid] = { getSpeciesName(t,nid) };
                //ti.enets[nid] = { getSpeciesName(t,nid) };


                for ( auto cid : t->getSonsId(nid) ) {
                    prepareTree( t, ti, cid );
                }

                for ( auto cid : t->getSonsId(nid) ) {
                    for ( auto l : ti.leaves[cid] ) { ti.leaves[nid].insert(l); }
                    for ( auto l : ti.subnodes[cid] ) { ti.subnodes[nid].insert(l); }
                    for ( auto l : ti.enets[cid] ) { ti.enets[nid].insert(l); }
                    for ( auto l : ti.subspec[cid] ) { ti.subspec[nid].insert(l); }
                }
            }
        }




    }

}

template bool Utils::appendNext<MultiOpt::dvsT, MultiOpt::QueueCmp<MultiOpt::dvsT>>(double, 
    const vector<size_t, StackAllocator<size_t>>&, 
    const vector<size_t, StackAllocator<size_t>>& , 
    vector<MultiOpt::dvsT>& ,
    MultiOpt::QueueCmp<MultiOpt::dvsT>& , 
    std::function< double(const vector<size_t, StackAllocator<size_t>>&) >& );


typedef boost::heap::skew_heap<MultiOpt::EdgeDerivation, boost::heap::compare<MultiOpt::CountedDerivCmp<MultiOpt::EdgeDerivation>>> heapT;
template bool Utils::appendNextWithEdge<heapT>( const size_t&,
                                                const vector<size_t, StackAllocator<size_t>>&,
                                                const vector<size_t, StackAllocator<size_t>>& ,
                                                heapT&,
                                                std::function< double(const size_t&, const vector<size_t, StackAllocator<size_t>>&) >&);
/*
template bool Utils::appendNextWithEdgeOrig<heapT>( const size_t&,
                                                const vector<size_t>&,
                                                const vector<size_t>& ,
                                                heapT&,
                                                std::function< double(const size_t&, const vector<size_t>&) >&);
*/
