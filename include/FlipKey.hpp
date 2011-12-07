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
    typedef tuple<int,int> NodeTupT;

public:
    friend ostream& operator<<(ostream& output, const FlipKey& flip) {
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

    FlipKey( int u, int v, bool f, bool r ) :
        _nodes( u < v ? NodeTupT(u,v) : NodeTupT(v,u)  ),
        _f(f), _r(r) {}

    FlipKey( const FlipKey& o ) :
        _nodes( NodeTupT(o.u(), o.v()) ), _f(o.f()), _r(o.r()) {}


    inline bool operator == (const FlipKey& other) const {
        return (_nodes == other._nodes) && (_f == other._f && _r == other._r);
    }

    inline bool operator != (const FlipKey& other) const {
        return !(*this == other);
    }

    int arity() const {
        if ( get<0>(_nodes) == get<1>(_nodes) ) { return 1; } else { return 2; }
    }

    const tuple<bool,bool> getDirTuple() const { return make_tuple(_f, _r); }

    int u() const { return get<0>(_nodes);}
    int v() const { return get<1>(_nodes);}
    bool f() const { return _f; }
    bool r() const { return _r; }

    std::size_t hashCode() const {
        std::hash<int> h;
        return (h( get<0>(_nodes) ) + h( get<1>(_nodes) )) + ( _f ? 2 : 0 ) + ( _r ? 1 : 0 );
    }

 private:
    NodeTupT _nodes;
    bool _f, _r;
};

    namespace std {
        template<>
        class hash<FlipKey> {
        public:
            std::size_t operator()(const FlipKey& k) const {
                return k.hashCode();
            }
        };
    }

    namespace boost {
        template<>
        class hash<FlipKey> {
        public:
            std::size_t operator()(const FlipKey& k) const {
                return k.hashCode();
            }
        };
    }


FlipKey flipBoth( const FlipKey& k ) { return FlipKey( k.u(), k.v(), !k.f(), !k.r() ); }
FlipKey flipForward( const FlipKey& k ) { return FlipKey( k.u(), k.v(), !k.f(), k.r() ); }
FlipKey flipReverse( const FlipKey& k ) { return FlipKey( k.u(), k.v(), k.f(), !k.r() ); }

#endif /** ARG_PARSER_HPP */
