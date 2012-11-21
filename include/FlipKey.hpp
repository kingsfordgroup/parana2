#ifndef FLIP_KEY_HPP
#define FLIP_KEY_HPP

#include <tuple>
#include <ostream>
#include <functional>

using std::ostream;
using std::tuple;
using std::make_tuple;
using std::string;
using std::get;


class FlipKey {
public:        
    typedef tuple<int,int> NodeTupT;
    typedef int NodeIndexT;

    friend ostream& operator<<(ostream& output, const FlipKey& flip);

    FlipKey( int u, int v, bool f, bool r );
    FlipKey( const FlipKey& o );

    bool operator == (const FlipKey& other) const;

    bool operator != (const FlipKey& other) const;

    int arity() const;

    const tuple<bool,bool> getDirTuple() const;

    int u() const;
    int v() const;
    bool f() const;
    bool r() const;

    std::size_t hashCode() const;

private:
    NodeTupT _nodes;
    bool _f, _r;
};


FlipKey flipBoth( const FlipKey& k );
FlipKey flipForward( const FlipKey& k );
FlipKey flipReverse( const FlipKey& k );

namespace boost {
    template<>
    class hash<FlipKey> {
    public:
        std::size_t operator()(const FlipKey& k) const;
    };
}

namespace std {
    template<>
    class hash<FlipKey> {
    public:
        size_t operator()(const FlipKey& k) const;
    };
}


#endif /** ARG_PARSER_HPP */
