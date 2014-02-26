#include <bitset>
#include <cassert>
#include <cstdint>
#include <iostream>
using namespace std;

const int N = 20;
const int M = N * N;

/*
  We precompute a lookup table for the game of life function for
  speed's sake.
 */
const int LIFE_TABLE = 1 << 9;
bool life_step[LIFE_TABLE];

void init( ) {
	for( int i = 0; i < LIFE_TABLE; i++ ) {
		int t = __builtin_popcount(i);
		if( t == 3 )
			life_step[i] = true;
		if( t == 4 )
			life_step[i] = (i & (1 << 4));
	}
}

/*
  We have a variety of helper tools to deal with a very fundamental
  data representation: a grid, especially 20x20 or 5x5.  Onto this
  grid we put a coordinate system, so that the center is (0,0).  We
  assume that this cell is "inside" the boundaries of our 20x20 grid,
  though this may be an assumption that will later change.

  The x- and y-coordinates range in a manner relative to the "center"
  cell of the grid, eg they could go from -2 to 2.  The orientation of
  the rows is so that y=-2 would represent the *top* row when drawn on
  paper.

  If the boundary of the grid is within the XxY cell, that information
  is represented.  Specifically, a positive dx indicates how many
  columns on the left are outside of the grid, while a negative the
  opposite on the other side.  A dx of 0 means at least part of each
  column is inside.  The same for dy, where positive dy means rows
  missing at the top.  Example:

  X....  X = outside
  X....  . = inside
  X....
  XXXXX  dx = +1  cx = 0
  XXXXX  dy = -2  cy = 0

  One relatively important thing is that every 5x5 such grid can be
  transformed into a 26-bit representation, which will be used to
  index into an array.
 */

template <int X, int Y>
struct grid {
	// (0,0) in our x-y coordinate system corresponds to (cx,cy)
	// in the usual system with (0,0) at the corner
	int cx, cy;

	// this represents how much of the grid is "out of bounds", in
	// the manner described above
	int dx, dy;

	bitset<X*Y> g;

	grid( ) : cx(0), cy(0), dx(0), dy(0), g(0) {
	}

	bool represented_uncentered( int x, int y ) const {
		return 0 <= x and x < X and 0 <= y and y < Y;
	}
	bool represented( int x, int y ) const {
		return represented_uncentered( x+cx, y+cy );
	}

	bool inside_uncentered( int x, int y ) const {
		return represented_uncentered(x,y) and
			0 <= (x-dx) and (x-dx) < X and 0 <= (y-dy) and (y-dy) < Y;
	}
	bool inside( int x, int y ) const {
		return inside_uncentered(x+cx, y+cy);
	}

	int get_int_uncentered( int x, int y ) const {
		if( not inside_uncentered(x,y) )
			return -1;
		return get_raw_uncentered(x,y);
	}
	int get_int( int x, int y ) const {
		return get_int_uncentered(x+cx, y+cy);
	}
	bool get_bool( int x, int y ) const {
		return ((not inside(x,y)) and get_raw(x,y));
	}
	bool get_raw_uncentered( int x, int y ) const {
		return g[ x*Y + y ];
	}
	bool get_raw( int x, int y ) const {
		return get_raw_uncentered( x+cx, y+cy );
	}

	void set_uncentered( int x, int y, bool v ) {
		g.set( x*Y + y, v );
	}
	void set( int x, int y, bool v ) {
		set_uncentered( x+cx, y+cy, v );
	}
};
typedef uint32_t encoding;

template <int X, int Y>
bool operator != ( const grid<X,Y>& g, const grid<X,Y>& h ) {
	if( g.dx != h.dx )
		return true;
	if( g.dy != h.dy )
		return true;
	return g.g != h.g;
}

/*
  The main encoding is fairly simple.  It does not encode (cx,cy).  So
  it just lists, in order, the different possible (dx,dy,grid) tuples,
  beginning with (dx,dy) = (0,0) and afterwards proceeding in
  lexicographic order -- eg, for a centered 5x5 grid, beginning with
  (-2,-2) and (-2,-1).  See the code for more details.
 */
template <int X, int Y>
encoding encode( grid<X,Y> g ) {
	// Optimized for dx = dy = 0
	if( g.dx == 0 and g.dy == 0 ) {
		return g.g.to_ulong();
	}

	// Otherwise, build up to number
	encoding e = (1 << (X*Y));
	for( int dx = -(X-1); dx <= (X-1); dx++ ) {
	for( int dy = -(Y-1); dy <= (Y-1); dy++ ) {
		if( dx == 0 and dy == 0 )
			continue;
		
		int width  = X - abs(dx);
		int height = Y - abs(dy);
		int area   = width * height;
		if( dx != g.dx or dy != g.dy ) {
			e += (1 << area);
		} else {
			encoding temp = 0;
			int k = 0;
			for( int x = 0; x < X; x++ ) {
			for( int y = 0; y < Y; y++ ) {
				int t = g.get_int_uncentered(x,y);
				if( t >= 0 ) {
					if( t == 1 )
						temp |= 1 << k;
					k++;
				}
			}
			}
			e += temp;
			return e;
		}
	}
	}
	// If we're here, then this better represent the empty grid
	assert( g.dx == X and g.dy == Y );
	return e;
}
template <int X, int Y>
grid<X,Y> decode( encoding f ) {
	grid<X,Y> g;
	// Optimized for dx = dy = 0
	if( f < (1 << (X*Y)) ) {
		g.dx = 0;
		g.dy = 0;
		g.g = bitset<X*Y>(f);
		return g;
	}
	
	// Otherwise, build up to number
	encoding e = (1 << (X*Y));
	for( int dx = -(X-1); dx <= (X-1); dx++ ) {
	for( int dy = -(Y-1); dy <= (Y-1); dy++ ) {
		if( dx == 0 and dy == 0 )
			continue;
		
		int width  = X - abs(dx);
		int height = Y - abs(dy);
		int area   = width * height;
		if( f >= e + (1 << area) ) {
			e += (1 << area);
		} else {
			g.dx = dx;
			g.dy = dy;

			encoding temp = f - e;
			int k = 0;
			for( int x = 0; x < X; x++ ) {
			for( int y = 0; y < Y; y++ ) {
				if( g.inside_uncentered(x,y) ) {
					g.set_uncentered( x,y, (temp & (1 << k)) != 0 );
					k++;
				}
			}
			}
			return g;
		}
	}
	}
	assert( e == f );
	g.dx = X;
	g.dy = Y;
	return g;
}

void test( int i ) {
	if( (i % 1000000) == 0 )
		cout << i/1000000 << " million" << endl;

	encoding  e = i;
	grid<5,5> g = decode<5,5>( e );
	encoding  f = encode<5,5>( g );
	grid<5,5> h = decode<5,5>( f );

	if( e != f ) {
		cout << e << " " << f << "\n";
		cout << "Failed on encoding!\n";
		exit(10);
	}
	if( g != h ) {
		cout << e << " " << f << "\n";
		cout << "Failed on grids!\n";
		exit(10);
	}
}

int main( ) {
	init();

	// Test encoding
	int i;
	for( i = 0; i < (1 << 25); i++ ) {
		test(i);
	}
	for( ; i < 38183849; i++ ) {
		test( i );
	}

	return 0;
}
