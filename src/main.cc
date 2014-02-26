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
  data representation: a 5x5 grid.  Onto this grid we put a coordinate
  system, so that the center is (0,0).  We assume that this cell is
  "inside" the boundaries of our 20x20 grid, though this may be an
  assumption that will later change.

  The x- and y-coordinates go from -2 to 2, with y=-2 representing the
  *top* row when drawn on paper.

  If the boundary of the grid is within the 5x5 cell, that information
  is represented.  Specifically, a positive dx indicates how many
  columns on the left are outside of the grid, while a negative the
  opposite on the other side.  A dx of 0 means at least part of each
  column is inside.  The same for dy, where positive dy means rows
  missing at the top.  Example:

  X....  X = outside
  X....  . = inside
  X....
  XXXXX  dx = +1
  XXXXX  dy = -2

  One relatively important thing is that every such grid can be
  transformed into a 26-bit representation, which will be used to
  index into an array.
 */
struct grid25 {
	int dx, dy;
	bitset<25> g;

	bool inside( int x, int y ) const {
		return -2 <= x and x <= 2 and -2 <= y and y <= 2 and 
			-2 <= (x-dx) and (x-dx) <= 2 and -2 <= (y-dy) and (y-dy) <= 2;
	}

	int get_int( int x, int y ) const {
		if( not inside(x,y) )
			return -1;
		return get_raw(x,y);
	}
	bool get_bool( int x, int y ) const {
		return ((not inside(x,y)) and get_raw(x,y));
	}
	bool get_raw( int x, int y ) const {
		return g[ 5*(x+2) + (y+2) ];
	}

	void set( int x, int y, bool v ) {
		g.set( 5*(x+2) + (y+2), v );
	}
};
typedef int32_t encoding26;

bool operator != ( const grid25& g, const grid25& h ) {
	if( g.dx != h.dx )
		return true;
	if( g.dy != h.dy )
		return true;
	for( int x = -2; x <= 2; x++ ) {
		for( int y = -2; y <= 2; y++ ) {
			if( g.get_bool(x,y) != h.get_bool(x,y) )
				return true;
		}
	}
	return false;
}

/*
  The encoding is 26 bits.  If the first (most significant bit) is 0,
  then the remaining 25 bits are just the grid in the lexicographic
  order beginning with (-2,-2) and (-2,-1).

  If the first bit is 1, then the next five bits encode (dx,dy) in
  lexicographic order from 0 meaning (-2,-2) and 1 meaning (-2,-1) to
  24 meaning (2,2).  It is not allowed for these bits to be 12 meaning
  (0,0) or be above 24.

  The remaining 20 bits represent the grid.  For the particular order
  of these bits, as well as the order of the five bits mentioned
  above, see the code.
 */
encoding26 encode( grid25 g ) {
	if( g.dx == 0 and g.dy == 0 ) {
		return g.g.to_ulong();
	}

	bitset<26> temp;
	temp[25] = true;
	// Encode upper five bits
	unsigned int dxy = 5 * (g.dx + 2) + (g.dy + 2);
	for( int i = 0; i < 5; i++ ) {
		temp.set( 20+i, dxy & (1 << i) );
	}

	// Encode rest of the bits
	int k = 0;
	for( int x = -2; x <= 2; x++ ) {
		for( int y = -2; y <= 2; y++ ) {
			int t = g.get_int(x,y);
			if( t >= 0 ) {
				temp.set( k, t );
				k++;
			}
		}
	}
	return temp.to_ulong();
}
grid25 decode( encoding26 e ) {
	grid25 g;
	if( (e & (1 << 25)) == 0 ) {
		g.dx = 0;
		g.dy = 0;
		g.g = bitset<25>(e);
		return g;
	}

	bitset<26> temp( e );
	// Read top five bits
	unsigned int dxy = 0;
	for( int i = 0; i < 5; i++ ) {
		if( temp[20+i] )
			dxy |= 1 << i;
	}
	assert( dxy < 25 );
	g.dx = (dxy / 5) - 2;
	g.dy = (dxy % 5) - 2;
	assert( g.dx != 0 or g.dy != 0 );

	// Read rest of the bits
	int k = 0;
	for( int x = -2; x <= 2; x++ ) {
		for( int y = -2; y <= 2; y++ ) {
			if( g.inside(x,y) ) {
				g.set(x,y, temp[k] );
				k++;
			}
		}
	}
	return g;	
}

void test( int i ) {
	if( (i % 1000000) == 0 )
		cout << i/1000000 << " million" << endl;

	encoding26 e = i;
	grid25 g = decode( e );
	encoding26 f = encode( g );
	grid25 h = decode( f );

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
	for( int dx = -2; dx <= 2; dx++ ) {
		for( int dy = -2; dy <= 2; dy++ ) {
			if( dx == 0 and dy == 0 )
				continue;
			unsigned int dxy = 5 * (dx + 2) + (dy + 2);
			int width  = 5 - abs(dx);
			int height = 5 - abs(dy);
			int area   = width * height;
			for( int j = 0; j < (1 << area); j++ ) {
				i = (1 << 25);
				i += dxy * (1 << 20);
				i += j;
				test( i );
			}
		}
	}

	return 0;
}
