#ifndef UTILS_HPP
#define UTILS_HPP

#include <unordered_map>
#include <unordered_set>

namespace Utils {
    using std::unordered_map;
    using std::unordered_set;
    using std::string;

    class TreeInfo {
    public:
        typedef unordered_map<int, unordered_set<string>> NodeMapT;
        typedef unordered_map<int, unordered_set<int>> NodeIndMapT;
        NodeIndMapT subnodes;
        NodeIndMapT leaves;
        NodeMapT enets;
        string name;

        TreeInfo( const string n ) : name(n), subnodes(NodeIndMapT()),
                                     leaves(NodeIndMapT()),
                                     enets(NodeMapT()) { cout << "creating tree info: " << name << "\n"; }

        bool inSubnodesOf( int u,  int v ) {
            return subnodes[u].find(v) != subnodes[u].end();
        }
        bool inSubnodesOf( int u,  int v ) const {
            return subnodes.find(u) != subnodes.cend() && subnodes.find(u)->second.find(v) != subnodes.find(u)->second.cend();
        }


    };


    bool differentExtantNetworks( const TreeInfo& tinfo, int u, int v ) {
        if ( u == v ) { return false; }
        for (auto uit = tinfo.enets.find(u)->second.cbegin(); uit != tinfo.enets.find(u)->second.cend(); ++uit) {
            if ( tinfo.enets.find(v)->second.find(*uit) != tinfo.enets.find(v)->second.cend() ) {
                return false;
            }
        }
        return true;
    }

}

#endif // UTILS_HPP
