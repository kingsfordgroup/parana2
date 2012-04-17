#include <tuple>
#include <ostream>
#include <functional>
#include <boost/functional/hash.hpp>
#include "FlipKey.hpp"

FlipKey::FlipKey( int u, int v, bool f, bool r ) :
    _nodes( u < v ? NodeTupT(u,v) : NodeTupT(v,u)  ),
    _f(f), _r(r) {}

FlipKey::FlipKey( const FlipKey& o ) :
    _nodes( NodeTupT(o.u(), o.v()) ), _f(o.f()), _r(o.r()) {}


bool FlipKey::operator == (const FlipKey& other) const {
    return (_nodes == other._nodes) && (_f == other._f && _r == other._r);
}

bool FlipKey::operator != (const FlipKey& other) const {
    return !(*this == other);
}

int FlipKey::arity() const {
    if ( get<0>(_nodes) == get<1>(_nodes) ) { return 1; } else { return 2; }
}

const tuple<bool,bool> FlipKey::getDirTuple() const { return make_tuple(_f, _r); }

int FlipKey::u() const { return get<0>(_nodes);}
int FlipKey::v() const { return get<1>(_nodes);}
bool FlipKey::f() const { return _f; }
bool FlipKey::r() const { return _r; }

std::size_t FlipKey::hashCode() const {
    std::hash<int> h;
    return (h( get<0>(_nodes) ) + h( get<1>(_nodes) )) + ( _f ? 2 : 0 ) + ( _r ? 1 : 0 );
}

FlipKey flipBoth( const FlipKey& k ) { return FlipKey( k.u(), k.v(), !k.f(), !k.r() ); }
FlipKey flipForward( const FlipKey& k ) { return FlipKey( k.u(), k.v(), !k.f(), k.r() ); }
FlipKey flipReverse( const FlipKey& k ) { return FlipKey( k.u(), k.v(), k.f(), !k.r() ); }

ostream& operator<<(ostream& output, const FlipKey& flip) {
        string flipStr;
        if (flip._f && flip._r) {
            flipStr = "<-->";
        } else if ( flip._f ) {
            flipStr = "-->";
        } else if ( flip._r ) {
            flipStr = "<--";
        } else {
            flipStr = "X";
        }

        output << "[" << get<0>(flip._nodes) << ", " << get<1>(flip._nodes) <<  "]" << " : " << flipStr;
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
