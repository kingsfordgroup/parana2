#include "TreeUtils.hpp"
#include "MultiOpt.hpp"

namespace Utils {


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
                     const vector<size_t>& inds,
                     const vector<size_t>& sizes,
                     vector<pqT>& pq,
                     pqCompT& pqComp,
                     std::function< double(const vector<size_t>&) >& computeScore ) {

        size_t i = 0;
        while ( i < inds.size() ) {
            vector<size_t> newInds(inds);
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

    typedef tuple<double, vector<size_t> > dvsT;
    class QueueCmp {
    public:
        bool operator() ( const dvsT& lhs, const dvsT& rhs ) {
            return get<0>(lhs) > get<0>(rhs);
        }
    };



    namespace Trees {
        typedef unique_ptr<bpp::TreeTemplate<bpp::Node>> TreePtrT;


        bool differentExtantNetworks( const TreeInfo& tinfo, int u, int v ) {
            if ( u == v ) { return false; }
            for (auto uit = tinfo.enets.find(u)->second.cbegin(); uit != tinfo.enets.find(u)->second.cend(); ++uit) {
                if ( tinfo.enets.find(v)->second.find(*uit) != tinfo.enets.find(v)->second.cend() ) {
                    return false;
                }
            }
            return true;
        }


        void labelTree( const TreePtrT& tree ) {
            for ( auto nid : tree->getNodesId() ) {
                for ( auto prop : tree->getBranchPropertyNames(nid) ) {
                    tree->setNodeName(nid, reinterpret_cast<bpp::BppString*>(tree->getBranchProperty(nid,prop))->toSTL() );
                }
            }
        }

        TreePtrT readNewickTree( const std::string& treeName ) {
            unique_ptr<bpp::Newick> newickReader( new bpp::Newick(true,true) ); //No comment allowed!
            newickReader->enableExtendedBootstrapProperty("name");
            try {
                TreePtrT tree( newickReader->read(treeName) ); // Tree in file
                cout << "Tree has " << tree->getNumberOfNodes() << " nodes." << endl;
                return tree;
            } catch(std::exception& e) {
                cerr << "Caught Exception: [" << e.what() << "]\n";
                abort();
            }
        }


        const std::string getName( TreePtrT& t, int nid) {
            if (t->hasNodeName(nid)) {
                return t->getNodeName(nid);
            } else {
                return reinterpret_cast<bpp::BppString*>( t->getNodeProperty(nid,"name"))->toSTL();
            }
        }

        string getExtantNetwork(const string& s) {
            return ( s.find("LOST") != string::npos ) ? "LOST" : s.substr( s.find_last_of('_') );
        };


        void prepareTree( TreePtrT& t, TreeInfo& ti, int nid ) {
            // if the current node is not the root
            if ( nid != t->getRootId() ) {
                auto fid = t->getFatherId(nid);
                auto parentDeathT = get<1>(ti.extantInterval[fid]);
                assert ( ti.extantInterval.find(fid) != ti.extantInterval.end() );
                double fdist = 0.0;
                if ( t->isLeaf(nid) ) {
                    fdist = std::numeric_limits<double>::infinity();
                } else if ( t->hasDistanceToFather(nid) ) {
                    fdist = t->getDistanceToFather(nid);
                }
                ti.extantInterval[ nid ] = make_tuple( parentDeathT, parentDeathT + fdist );
            }


            if (t->isLeaf(nid)) {
                ti.leaves[nid] = {nid};
                ti.subnodes[nid] = {nid};
                string enet = getExtantNetwork(getName(t,nid));
                ti.enets[nid] = { };
                if( enet != "LOST" ) { ti.enets[nid].insert(enet); }
            } else {
                ti.leaves[nid] = { nid };
                ti.subnodes[nid] = { nid };
                ti.enets[nid] = unordered_set<string>();


                for ( auto cid : t->getSonsId(nid) ) {
                    prepareTree( t, ti, cid );
                }

                for ( auto cid : t->getSonsId(nid) ) {
                    for ( auto l : ti.leaves[cid] ) { ti.leaves[nid].insert(l); }
                    for ( auto l : ti.subnodes[cid] ) { ti.subnodes[nid].insert(l); }
                    for ( auto l : ti.enets[cid] ) { ti.enets[nid].insert(l); }
                }
            }
        }




    }

}

template bool Utils::appendNext<MultiOpt::dvsT, MultiOpt::QueueCmp>( double, const vector<size_t>&, const vector<size_t>& , vector<MultiOpt::dvsT>& , MultiOpt::QueueCmp& , std::function< double(const vector<size_t>&) >& );
