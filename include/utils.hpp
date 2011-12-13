#ifndef UTILS_HPP
#define UTILS_HPP

#include <unordered_map>
#include <unordered_set>
#include <tuple>

namespace Utils {
    using std::unordered_map;
    using std::unordered_set;
    using std::string;
    using std::tuple;

    class TreeInfo {
    public:
        typedef unordered_map<int, unordered_set<string>> NodeMapT;
        typedef unordered_map<int, unordered_set<int>> NodeIndMapT;
        typedef unordered_map<int, tuple<double,double>> NodeIntervalMapT;
        NodeIndMapT subnodes;
        NodeIndMapT leaves;
        NodeMapT enets;
        NodeIntervalMapT extantInterval;
        string name;

        TreeInfo( const string n ) : name(n), subnodes(NodeIndMapT()),
                                     leaves(NodeIndMapT()),
                                     enets(NodeMapT()),
                                     extantInterval(NodeIntervalMapT()) { cout << "creating tree info: " << name << "\n"; }

        bool inSubnodesOf( int u,  int v ) {
            return subnodes[u].find(v) != subnodes[u].end();
        }
        bool inSubnodesOf( int u,  int v ) const {
            return subnodes.find(u) != subnodes.cend() && subnodes.find(u)->second.find(v) != subnodes.find(u)->second.cend();
        }

        double intervalDistance( int u, int v ) {
            //return static_cast<const TreeInfo&>(*this).intervalDistance(u,v);
            auto birthU = get<0>( extantInterval[u] );
            auto deathU = get<1>( extantInterval[u] );
            auto birthV = get<0>( extantInterval[v] );
            auto deathV = get<1>( extantInterval[v] );
            // If the ranges overlap, the distance is 0
            if ( (birthU <= deathV) && (deathU >= birthV) ) { return 0.0; }
            // otherwise
            if ( birthU >= deathV ) {
                return birthU - deathV;
            }

            if ( birthV >= deathU ) {
                return birthV - deathU;
            }
        }
        /*
        double intervalDistance( int u, int v ) {
            return static_cast<const TreeInfo&>(*this).intervalDistance(u,v);
            }*/
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
