#ifndef MODEL_HPP
#define MODEL_HPP

#include <array>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <cmath>

#include <Bpp/Phyl/TreeTools.h>

#include "MultiOpt.hpp"
#include "TreeUtils.hpp"
#include "FlipKey.hpp"

class Model {
private:
	/**
	* Typedefs
	*/
	typedef double ProbabilityT;
	typedef std::tuple<bool,bool> FlipTupleT;
	//typedef std::unordered_map< FlipTupleT, std::unordered_map< FlipTupleT, ProbabilityT > > TransitionMapT;
	typedef std::array< std::array<ProbabilityT,4>, 4> TransitionMapT;

	const Utils::TreeInfo& _tinfo;
	Utils::Trees::TreePtrT& _t;
	TransitionMapT _prob;
	//std::unique_ptr<bpp::DistanceMatrix> _distMat;

	ProbabilityT _transitionProbability( const FlipKey& child, const FlipKey& parent, double dist ){
		// ignore distance right now
		// 0 + (0.92-0)/((1+0.25*e^{-5.0*{x-1.5}})^{1/1.})
		//auto u = _prob[ FlipTupleT(parent.f(),parent.r()) ][ FlipTupleT(child.f(), child.r())  ];
		bool keepingEdge = parent.f() and child.f();
		bool loosingEdge = parent.f() and (!child.f());
		bool gainingEdge = (!parent.f()) and child.f();
		bool keepingNoEdge  = (!parent.f()) and (!child.f());
		if ( loosingEdge ) {
			auto u = _prob[ FlipState::both ][ FlipState::none];
			return (u) / ((1.0 + 0.25 * std::exp( -5.0 * (dist - 1.5) ) ));
		} else if ( keepingEdge ) {
			auto u = _prob[ FlipState::both ][ FlipState::none];
			return 1.0 - ((u) / ((1.0 + 0.25 * std::exp( -5.0 * (dist - 1.5) ) )));
		} else if ( gainingEdge ) {
			auto u = _prob[ FlipState::none ][ FlipState::both];
			return (u)/((1.0+ std::exp(-3.0*(dist-1.5))));
		} else if ( keepingNoEdge ) {
			auto u = _prob[ FlipState::none ][ FlipState::both];
			return 1.0 - (u)/((1.0+ std::exp(-3.0*(dist-1.5))));
		}
	}
	
	void _readUndirectedModel( std::ifstream& ifile ) {
		// 0->0, 0->1, 1->0, 1->1
		std::vector<FlipState> tf{FlipState::none, FlipState::both};
		for ( const auto& e0 : tf ) {
			for ( const auto& e1 : tf ) {
				ifile >> _prob[e0][e1];
			}
		}
	}

	void _readDirectedModel( std::ifstream& ifile ) {
		// 00->00, 00->01, 00->10, 00->11
		// 01->00, 01->01, 01->10, 01->11
		// 10->00, 10->01, 10->10, 10->11
		// 11->00, 11->01, 11->10, 11->11
		std::vector<FlipState> tf{FlipState::none, FlipState::both};
		for ( const auto& f0 : tf ) {
			for ( const auto& r0 : tf ) {
				for ( const auto& f1 : tf ) {
					for ( const auto& r1 : tf ) {
						//ifile >> _prob[FlipTupleT(f0,r0)][FlipTupleT(f1,r1)];
						ifile >> _prob[f0][f1];
					}
				}
			}
		}
	}

	void _readModelFromFile( const std::string& mfile ) {
		using std::string;
		string line;
		std::ifstream ifile(mfile);		
		ifile >> line;
		std::cerr << "line = [" << line << "] \n";
		if ( line == "undirected" ) { _readUndirectedModel( ifile ); }
		else if ( line == "directed" ) { _readDirectedModel( ifile ); }
		else { std::cerr << "Model must be one of {undirected|directed}\n"; std::abort();}
		ifile.close();
	}

public:

	Model( const std::string& mfile, const Utils::TreeInfo& tinfo, Utils::Trees::TreePtrT& t ) : 
	    _tinfo(tinfo), _t(t) {
		std::cerr << "Reading from file\n";
		_readModelFromFile(mfile);
		//std::cerr << "computing distance matrix . . .";
		//_distMat.reset( bpp::TreeTools::getDistanceMatrix( *_t.get() ) );
		//std::cerr << "done\n";
	}
	
	ProbabilityT leafProbability( const FlipKey& leaf, const FlipKey& obs ) {
		return _transitionProbability( leaf, obs, 0.0 );
	}

	ProbabilityT transitionProbability( const FlipKey& child, const FlipKey& parent ) {
		FlipKey::NodeIndexT parentU, parentV, childU, childV;
		std::tie(parentU, parentV) = std::make_tuple(parent.u(), parent.v());
		std::tie(childU, childV) = std::make_tuple(child.u(), child.v());
	    
		FlipKey::NodeIndexT pnode, cnode, onode;
		if ( parentU == childU ) {
			pnode = parentU; cnode = childV; onode = parentV;
		} else if ( parentU == childV ) {
			pnode = parentU; cnode = childU; onode = parentV;
		} else if ( parentV == childU ) {
			pnode = parentV; cnode = childV; onode = parentU;
		} else if ( parentV == childV ) {
			pnode = parentV; cnode = childU; onode = parentU;
		}


		//auto dist = _t->getDistanceToFather(cnode);//1.0;
		//auto dist = bpp::TreeTools::getDistanceBetweenAnyTwoNodes( *_t.get(), onode, cnode);
		//auto dist = _distMat->operator()(onode, cnode);
		auto dist = 1.0;

		//auto dist = _tinfo.intervalDistance(pnode, cnode);

		return _transitionProbability( child, parent, dist );
	}
};

#endif // MODEL_HPP