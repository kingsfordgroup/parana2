#ifndef MODEL_HPP
#define MODEL_HPP

#include <array>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <cmath>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

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
	typedef std::array< ProbabilityT, 4 > PriorProbabilityMapT;

	const Utils::TreeInfo& _tinfo;
	Utils::Trees::TreePtrT& _t;
	TransitionMapT _dupProb;
	TransitionMapT _specProb;
	PriorProbabilityMapT _priorProb;
	//std::unique_ptr<bpp::DistanceMatrix> _distMat;

	ProbabilityT _transitionProbability( const FlipKey& child, const FlipKey& parent, double dist ){
		// ignore distance right now
		// 0 + (0.92-0)/((1+0.25*e^{-5.0*{x-1.5}})^{1/1.})
		//auto r = _prob[ parent.state() ][ child.state() ];
		//return r;

		auto ps = parent.state();
		auto cs = child.state();

		
		bool keepingEdge = (ps == FlipState::both) and (cs == FlipState::both);
		bool loosingEdge = (ps == FlipState::both) and (cs == FlipState::none);
		bool gainingEdge = (ps == FlipState::none) and (cs == FlipState::both);
		bool keepingNoEdge  = (ps == FlipState::none) and (cs == FlipState::none);

		auto& _prob = (Utils::Trees::isSpeciationEvent(parent, _tinfo)) ? _specProb : _dupProb;
		return _prob[ps][cs];
		/*
		auto u = _prob[ps][cs];

		auto probLose = (u) / ((1.0 + 0.25 * std::exp( -5.0 * (dist - 1.5) ) ));
		auto probGain = (u) / ((1.0+ std::exp(-3.0*(dist-2.5))));

		if ( loosingEdge ) {
			return probLose;
		} else if ( keepingEdge ) {
			return 1.0 - probLose;
		} else if ( gainingEdge ) {
			return probGain;
		} else if ( keepingNoEdge ) {
			return 1.0 - probGain;
		}
		*/
	}
	
	void _readUndirectedModel( boost::property_tree::ptree& pt ) {
		
		auto none = FlipState::none;
		auto both = FlipState::both;		

		_dupProb[both][both] = pt.get<double>("model.dupProbs.keepInteraction");
		_dupProb[both][none] = 1.0 - _dupProb[both][both];
		_dupProb[none][both] = pt.get<double>("model.dupProbs.gainInteraction");
		_dupProb[none][none] = 1.0 - _dupProb[none][both];

		_specProb[both][both] = pt.get<double>("model.specProbs.keepInteraction");
		_specProb[both][none] = 1.0 - _specProb[both][both];
		_specProb[none][both] = pt.get<double>("model.specProbs.gainInteraction");
		_specProb[none][none] = 1.0 - _specProb[none][both];

		_priorProb[both] = pt.get<double>("model.priorProbs.ancestralEdgeExistence");
		_priorProb[none] = 1.0 - _priorProb[both];

		/*
		// 0->0, 0->1, 1->0, 1->1
		std::vector<FlipState> tf{FlipState::none, FlipState::both};
		for ( const auto& e0 : tf ) {
			for ( const auto& e1 : tf ) {
				ifile >> _prob[e0][e1];
			}
		}
		*/
	}

	void _readDirectedModel( boost::property_tree::ptree& pt ) {
		/*
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
		*/
	}

	void _readModelFromFile( const std::string& mfile ) {

		using std::string;
		using boost::property_tree::ptree;
		ptree pt;

		read_xml(mfile, pt);
		if ( pt.get<std::string>("model.type") == "undirected" ) {
			_readUndirectedModel(pt);
		} else if ( pt.get<std::string>("model.type") == "directed" ) {
			_readDirectedModel(pt);
		} else {
			std::cerr << "model type must be one of {undirected, directed}\n";
			std::exit(1);
		}

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
	
	inline ProbabilityT priorProbability(FlipState s) {
		return _priorProb[s];
	}

	ProbabilityT leafProbability( const FlipKey& leaf, const FlipKey& obs ) {
		return _transitionProbability( leaf, obs, 0.0 );
	}

	ProbabilityT transitionProbability( const FlipKey& child, const FlipKey& parent ) {
		FlipKey::NodeIndexT parentU, parentV, childU, childV;

		std::tie(parentU, parentV) = std::make_tuple(parent.u(), parent.v());
		std::tie(childU, childV) = std::make_tuple(child.u(), child.v());

		FlipKey::NodeIndexT pnode, cnode, onode;

		if (parentU == parentV) {
			auto dist =  _t->getDistanceToFather(childU) + _t->getDistanceToFather(childV);
			return _transitionProbability(child, parent, dist);
		} else {

			if ( parentU == childU ) {
				pnode = parentU; cnode = childV; onode = parentV;
			} else if ( parentU == childV ) {
				pnode = parentU; cnode = childU; onode = parentV;
			} else if ( parentV == childU ) {
				pnode = parentV; cnode = childV; onode = parentU;
			} else if ( parentV == childV ) {
				pnode = parentV; cnode = childU; onode = parentU;
			}

			auto dist = 1.0;
  		//auto dist = _t->getDistanceToFather(cnode);//1.0;
  		//auto dist = bpp::TreeTools::getDistanceBetweenAnyTwoNodes( *_t.get(), onode, cnode);
	  	//auto dist = _distMat->operator()(onode, cnode);
		  //auto dist = 1.0;
		  //auto dist = _tinfo.intervalDistance(cnode, pnode);

  		return _transitionProbability( child, parent, dist );
  	}

	}
};

#endif // MODEL_HPP