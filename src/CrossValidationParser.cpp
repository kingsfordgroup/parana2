//#include "CrossValidationTest.hpp"
#include "CrossValidationParser.hpp"
#include <deque>
#include <tuple>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <pugixml.hpp>

CrossValidationTestParser::CrossValidationTestParser( const std::string& fname ) : fname_(fname) {}

CrossValidationTest CrossValidationTestParser::parse() {
	pugi::xml_document doc;
	pugi::xml_parse_result result = doc.load_file( fname_.c_str() );
	std::cerr << "Reading cross validation sets from : " << fname_ << "\n";

	using std::make_tuple;
	using std::tie;
	
	CrossValidationTest crossValTests;

	if (result){
		auto cvtest = doc.child("cvtest");
		crossValTests.testSetName = cvtest.attribute("name").value();
		std::cerr << "testSetName = " << crossValTests.testSetName << "\n";

		// Add all of the cross-validation sets for this test
		for (pugi::xml_node cvset = cvtest.child("testset"); cvset; cvset = cvset.next_sibling("testset")) {
			crossValTests.cvSets.emplace_back( cvset.attribute("name").value() );

			// Add all the edges for this cross-validation set
			for (pugi::xml_node edge = cvset.child("edge"); edge; edge = edge.next_sibling("edge")) {
				crossValTests.cvSets.back().edges.emplace_back( 
				              edge.attribute("u").value(), edge.attribute("v").value() );
		    }
		}

	}

	return crossValTests;
}

/*
int main(int argc, char* argv[] ) {
	CrossValidationTestParser parser( argv[1] );
	parser.parse();
}
*/