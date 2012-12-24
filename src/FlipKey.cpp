#include <tuple>
#include <ostream>
#include <iostream>
#include <functional>
#include <boost/functional/hash.hpp>
#include "FlipKey.hpp"

FlipKey::FlipKey( int u, int v, FlipState s, bool connectU, bool connectV ) :
    _nodes( u < v ? NodeTupT(u,v) : NodeTupT(v,u)  ),
    _state(s), _connectU(connectU), _connectV(connectV) {
        if ( u > v ) { std::cerr << "swapping connect states\n"; std::swap(_connectU, _connectV); }
    }

FlipKey::FlipKey( int u, int v, bool f, bool r, bool connectU, bool connectV ) :
    _nodes( u < v ? NodeTupT(u,v) : NodeTupT(v,u)  ),
    _connectU(connectU), _connectV(connectV) {
        if ( f and r ){ _state = FlipState::both; }
        if ( f and !r ){ _state = FlipState::forward; }
        if ( !f and r ){ _state = FlipState::reverse; }
        if ( !f and !r ){ _state = FlipState::none; }
        if ( u > v ) { std::cerr << "swapping connect states\n"; std::swap(_connectU, _connectV); }
    }

FlipKey::FlipKey( const FlipKey& o ) :
    _nodes( NodeTupT(o.u(), o.v()) ), _state(o.state()), _connectU(o.connectU()), _connectV(o.connectV()) {}


bool FlipKey::operator == (const FlipKey& other) const {
    return (_nodes == other._nodes) && (_state == other._state) && 
           (_connectU == other._connectU) && (_connectV == other._connectV);
}

bool FlipKey::operator != (const FlipKey& other) const {
    return !(*this == other);
}

int FlipKey::arity() const {
    if ( get<0>(_nodes) == get<1>(_nodes) ) { return 1; } else { return 2; }
}

//const tuple<bool,bool> FlipKey::getDirTuple() const { return make_tuple(_f, _r); }

int FlipKey::u() const { return get<0>(_nodes);}
int FlipKey::v() const { return get<1>(_nodes);}
bool FlipKey::f() const { return (_state == FlipState::forward or _state == FlipState::both); }
bool FlipKey::r() const { return (_state == FlipState::reverse or _state == FlipState::both); }
FlipState FlipKey::state() const { return _state; }
bool FlipKey::connectU() const { return _connectU; }
bool FlipKey::connectV() const { return _connectV; }

std::size_t FlipKey::hashCode() const {
    std::hash<int> h;
    return (h( get<0>(_nodes) ) + h( get<1>(_nodes) )) + _state  + 
           (_connectU ? 2 : 0) + (_connectV ? 1 : 0);
}

FlipKey flipBoth( const FlipKey& k ) { return FlipKey( k.u(), k.v(), fsBoth[k.state()], k.connectU(), k.connectV() ); }
FlipKey flipForward( const FlipKey& k ) { return FlipKey( k.u(), k.v(), fsForward[k.state()], k.connectU(), k.connectV() ); }
FlipKey flipReverse( const FlipKey& k ) { return FlipKey( k.u(), k.v(), fsReverse[k.state()], k.connectU(), k.connectV() ); }

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

        output << "[" << get<0>(flip._nodes) << ", " << get<1>(flip._nodes) <<  "]" << " : " << flipStr << 
                  "(" << flip._connectU << ", " << flip._connectV << ")";
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
