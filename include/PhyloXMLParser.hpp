#include "TreeUtils.hpp"

class PhyloXMLParser {
public:
	PhyloXMLParser( const std::string& fname );
	Utils::Trees::TreePtrT::element_type* parse();
private:
	std::string _fname;
};