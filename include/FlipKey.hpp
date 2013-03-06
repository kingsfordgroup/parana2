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

enum class MustConnect: uint32_t {
    none=0,
    left=1,
    right=2,
    both=3
};

enum FlipState: uint32_t {
    none=0,
    forward=1,
    reverse=2,
    both=3
};

static const FlipState fsBoth[] = { 
    FlipState::both, FlipState::reverse, FlipState::forward, FlipState::none };

static const FlipState fsForward[] = { 
    FlipState::forward, // X  : ->
    FlipState::none,    // -> : X
    FlipState::both,    // <- : <->
    FlipState::reverse  // <->: <-
};

static const FlipState fsReverse[] = { 
    FlipState::reverse, // X  : <-
    FlipState::both,    // -> : <->
    FlipState::none,    // <- : X
    FlipState::forward  // <->: ->
};

void swapEndpoints( MustConnect& c );

class FlipKey {
public:        
    typedef tuple<int,int> NodeTupT;
    typedef int NodeIndexT;

    friend ostream& operator<<(ostream& output, const FlipKey& flip);

    FlipKey( int u, int v, FlipState t, MustConnect c=MustConnect::none );
    FlipKey( int u, int v, FlipState t, bool connectU, bool connectV );
    FlipKey( int u, int v, bool f, bool r, bool connectU=false, bool connectV=false );
    FlipKey( const FlipKey& o );

    bool operator == (const FlipKey& other) const;

    bool operator != (const FlipKey& other) const;

    int arity() const;
    bool isSelfLoop() const;
    //const tuple<bool,bool> getDirTuple() const;

    int u() const;
    int v() const;
    bool f() const;
    bool r() const;
    FlipState state() const;
    MustConnect mustConnect() const;
    bool connectU() const;
    bool connectV() const;

    std::size_t hashCode() const;

private:
    NodeTupT _nodes;
    FlipState _state;
    //bool _f, _r;
    //bool _connectU, _connectV;
    MustConnect _connect;
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
