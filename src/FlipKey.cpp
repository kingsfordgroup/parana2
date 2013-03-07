#include <tuple>
#include <ostream>
#include <iostream>
#include <functional>
#include <boost/functional/hash.hpp>
#include "FlipKey.hpp"

void swapEndpoints( MustConnect& c ) {
    switch (c) {
        case MustConnect::left :
          c = MustConnect::right;
        case MustConnect::right :
          c = MustConnect::left;
    }
}

FlipKey::FlipKey( int u, int v, FlipState s, MustConnect c ) :
    _nodes( u < v ? NodeTupT(u,v) : NodeTupT(v,u)  ),
    _state(s), _connect(c) {
        if ( u > v ) { std::cerr << "swapping connect states\n"; swapEndpoints(_connect); }
    }

FlipKey::FlipKey( int u, int v, FlipState s, bool connectU, bool connectV ) :
    _nodes( u < v ? NodeTupT(u,v) : NodeTupT(v,u)  ),
    _state(s) {
        if ( u > v ) { std::cerr << "swapping connect states\n"; std::swap(connectU, connectV); }
        if ( connectU and connectV ) { _connect = MustConnect::both; }
        else if ( connectU and !connectV ) { _connect = MustConnect::left; }
        else if ( !connectU and connectV ) { _connect = MustConnect::right; }
        else if ( !connectU and !connectV ) { _connect = MustConnect::none; }
    }

FlipKey::FlipKey( int u, int v, bool f, bool r, bool connectU, bool connectV ) :
    _nodes( u < v ? NodeTupT(u,v) : NodeTupT(v,u)  ) {
        if ( f and r ){ _state = FlipState::both; }
        else if ( f and !r ){ _state = FlipState::forward; }
        else if ( !f and r ){ _state = FlipState::reverse; }
        else if ( !f and !r ){ _state = FlipState::none; }
        
        if ( u > v ) { std::cerr << "swapping connect states\n"; std::swap(connectU, connectV); }
        if ( connectU and connectV ) { _connect = MustConnect::both; }
        else if ( connectU and !connectV ) { _connect = MustConnect::left; }
        else if ( !connectU and connectV ) { _connect = MustConnect::right; }
        else if ( !connectU and !connectV ) { _connect = MustConnect::none; }
    }

FlipKey::FlipKey( const FlipKey& o ) :
    _nodes( NodeTupT(o.u(), o.v()) ), _state(o.state()), _connect(o.mustConnect()) {}


bool FlipKey::operator == (const FlipKey& other) const {
    return (_nodes == other._nodes) && (_state == other._state) && (_connect == other._connect);
}

bool FlipKey::operator != (const FlipKey& other) const {
    return !(*this == other);
}

int FlipKey::arity() const {
    if ( get<0>(_nodes) == get<1>(_nodes) ) { return 1; } else { return 2; }
}

bool FlipKey::isSelfLoop() const {
    return arity()==1;
}

//const tuple<bool,bool> FlipKey::getDirTuple() const { return make_tuple(_f, _r); }

int FlipKey::u() const { return get<0>(_nodes);}
int FlipKey::v() const { return get<1>(_nodes);}
bool FlipKey::f() const { return (_state == FlipState::forward or _state == FlipState::both); }
bool FlipKey::r() const { return (_state == FlipState::reverse or _state == FlipState::both); }
FlipState FlipKey::state() const { return _state; }
MustConnect FlipKey::mustConnect() const { return _connect; }
bool FlipKey::connectU() const { 
    return (_connect == MustConnect::left or _connect == MustConnect::both); }

bool FlipKey::connectV() const { 
    return (_connect == MustConnect::right or _connect == MustConnect::both); }

std::size_t FlipKey::hashCode() const {
    size_t seed = 0;
    boost::hash_combine(seed, get<0>(_nodes));
    boost::hash_combine(seed, get<1>(_nodes));
    boost::hash_combine(seed, static_cast<int>(_state));
    boost::hash_combine(seed, static_cast<int>(_connect));
    return seed;
    //std::hash_combine<int> h;
    //return (h( get<0>(_nodes) ) + h( get<1>(_nodes) )) + static_cst<int>(_state)  + static_cast<int>(_connect);
}

FlipKey flipBoth( const FlipKey& k ) { return FlipKey( k.u(), k.v(), fsBoth[k.state()], k.mustConnect() ); }
FlipKey flipForward( const FlipKey& k ) { return FlipKey( k.u(), k.v(), fsForward[k.state()], k.mustConnect() ); }
FlipKey flipReverse( const FlipKey& k ) { return FlipKey( k.u(), k.v(), fsReverse[k.state()], k.mustConnect() ); }

ostream& operator<<(ostream& output, const FlipKey& flip) {
        string flipStr;
        switch ( flip._state ) {
            case FlipState::both:
              flipStr = "<-->"; break;
            case FlipState::forward:
              flipStr = "-->"; break;
            case FlipState::reverse:
              flipStr = "<--"; break;
            case FlipState::none:
              flipStr = "X"; break;
        }

        string connectStr;
        switch ( flip._connect ) {
            case MustConnect::both:
              connectStr = "(1,1)"; break;
            case MustConnect::left:
              connectStr = "(1,0)"; break;
            case MustConnect::right:
              connectStr = "(0,1)"; break;
            case MustConnect::none:
              connectStr = "(0,0)"; break;
        }

        output << "[" << get<0>(flip._nodes) << ", " << get<1>(flip._nodes) <<  "]" << " : " << flipStr << 
                  connectStr;
        return output;
    }

namespace boost {
    std::size_t hash<FlipKey>::operator()(const FlipKey& k) const {
        return k.hashCode();
    }
}

namespace std {
    size_t hash<FlipKey>::operator()(const FlipKey& k) const {
        return k.hashCode();
    }
}
