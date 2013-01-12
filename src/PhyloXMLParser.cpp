#include <deque>
#include <tuple>
#include <cstdlib>
#include <exception>

#include <pugixml.hpp>

#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Node.h>
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Exceptions.h>

#include "TreeUtils.hpp"
#include "PhyloXMLParser.hpp"

using Utils::Trees::TreePtrT;

PhyloXMLParser::PhyloXMLParser( const std::string& fname ) : _fname(fname) {}
	
TreePtrT::element_type* PhyloXMLParser::parse() {
	pugi::xml_document doc;
	pugi::xml_parse_result result = doc.load_file( _fname.c_str() );
	std::cerr << _fname << "\n";

	using std::make_tuple;
	using std::tie;
		
	auto tree = new TreePtrT::element_type;
	if (result){
		auto phyloxml = doc.child("phyloxml");
		auto phylogeny = phyloxml.child("phylogeny");

		auto root = phylogeny.child("clade");
		std::deque< std::tuple<pugi::xml_node, bpp::Node*> > toProcess;
		auto rootNode = new bpp::Node(root.child("name").child_value());
		tree->setRootNode(rootNode);
		toProcess.push_front( make_tuple(root, rootNode) );
		while ( !toProcess.empty() ) {
			pugi::xml_node clade;
			bpp::Node* node;
			std::tie(clade, node) = toProcess.front();
			toProcess.pop_front();


			auto nodeName = clade.child("name").child_value();
			//std::cerr << "NAME: " << nodeName << ", ";

			auto nodeDist = std::atof(clade.child("branch_length").child_value());
			auto species = clade.child("taxonomy").child("scientific_name").child_value();
			auto dup = clade.child("events").child("duplication") ? "T" : "F";

			//std::cerr << "dist: " << nodeDist << ", species: " << species << "\n";
			node->setDistanceToFather(nodeDist);

			bpp::BppString specName(species);
			node->setNodeProperty("S", specName);

			bpp::BppString dupEvent(dup);
			node->setNodeProperty("D", dupEvent);

			bpp::BppString nodeString(nodeName);
			node->setNodeProperty("GN", nodeString);

			//std::cerr << "children: [";
			for (pugi::xml_node subclade: clade.children("clade")) {

				auto nodeName = subclade.child("name").child_value();

				bpp::Node* subnode = new bpp::Node(nodeName);
				//std::cerr << nodeName << ",";

				node->addSon(subnode);
				toProcess.push_back( make_tuple(subclade, subnode) );
			}
			//std::cout << "]" << std::endl;
		}

		tree->resetNodesId();
		//simple_walker walker;
		//doc.traverse(walker);
	} else {
		std::cerr << "fail\n";
	}
	return tree;
}
