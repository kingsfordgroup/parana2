#ifndef __CROSSVALIDATIONPARSER__HPP__
#define __CROSSVALIDATIONPARSER__HPP__

#include <string>
#include <vector>

class StringEdge{
public:
	std::string u;
	std::string v;
	StringEdge( const std::string& _u, const std::string& _v ) : u(_u), v(_v) {}
};

class CrossValidationSet{
public:
	std::string name;
	std::vector< StringEdge > edges;

	CrossValidationSet( const std::string& _name ) : name(_name), edges(std::vector<StringEdge>()) {}
	CrossValidationSet( const std::string& _name, std::vector<StringEdge>& _edges ) : 
	                    name(_name), edges(std::move(_edges)) {}
};

class CrossValidationTest{
public:
	std::string testSetName;
	std::vector< CrossValidationSet > cvSets;
};

class CrossValidationTestParser {
public:
	CrossValidationTestParser( const std::string& fname );
	CrossValidationTest parse();
private:
	std::string fname_;
};

#endif // __CROSSVALIDATIONPARSER__HPP__