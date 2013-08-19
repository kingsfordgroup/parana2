#ifndef PARSIMONY_COSTS_HPP
#define PARSIMONY_COSTS_HPP

#include <array>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <cmath>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "MultiOpt.hpp"
#include "TreeUtils.hpp"
#include "FlipKey.hpp"


class ParsimonyCosts {
  public:
    using costRepT = std::tuple<double, std::string>;
    using costMapT = std::array<std::array<costRepT,4>,4>;
    using selfCostMapT = std::array<std::array<costRepT,2>,4>;

    ParsimonyCosts(const std::string& fname) : _leafCostRatio(1.0), _directed(false) { 
        _readCostModelFromFile(fname); 
    }

    void _readDirectedModel(boost::property_tree::ptree& pt) {}

    void _readUndirectedModel(boost::property_tree::ptree& pt) {
        auto none = FlipState::none;
        auto both = FlipState::both;    

        _dupCost[both][both] = make_tuple(pt.get<double>("model.dupCosts.keepInteraction"), "n");
        _dupCost[both][none] = make_tuple(pt.get<double>("model.dupCosts.looseInteraction"), "b-");
        _dupCost[none][both] = make_tuple(pt.get<double>("model.dupCosts.gainInteraction"), "b+");
        _dupCost[none][none] = make_tuple(pt.get<double>("model.dupCosts.keepNoInteraction"), "n");

        _specCost[both][both] = make_tuple(pt.get<double>("model.specCosts.keepInteraction"), "n");
        _specCost[both][none] = make_tuple(pt.get<double>("model.specCosts.looseInteraction"), "b-");
        _specCost[none][both] = make_tuple(pt.get<double>("model.specCosts.gainInteraction"), "b+");
        _specCost[none][none] = make_tuple(pt.get<double>("model.specCosts.keepNoInteraction"), "n");

    }    

    void _readCostModelFromFile(const std::string& fname) {
        using std::string;
        using boost::property_tree::ptree;
        ptree pt;

        read_xml(fname, pt);
        _leafCostRatio = pt.get<double>("model.leafCostRatio");
        
        if ( pt.get<std::string>("model.type") == "undirected" ) {
          _readUndirectedModel(pt);
          _directed = false;
        } else if ( pt.get<std::string>("model.type") == "directed" ) {
          _readDirectedModel(pt);
          _directed = true;
        } else {
          std::cerr << "cost model type must be one of {undirected, directed}\n";
          std::exit(1);
        }
    }

    double leafCostRatio() const { return _leafCostRatio; }
    bool isDirected() const { return _directed; }
    const costMapT& getCostDict() const { return _dupCost; }
    const costMapT& getSpecCostDict() const { return _specCost; }    

  private:
    costMapT _dupCost;
    costMapT _specCost;
    double _leafCostRatio;
    bool _directed;
};


#endif // PARSIMONY_COSTS_HPP