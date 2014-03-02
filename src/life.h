#ifndef LIFE_H
#define LIFE_H

#include "grid.h"

/*
  We precompute a lookup table for the game of life function for
  speed's sake.
 */
const int LIFE_TABLE = 1 << 9;
bool life_step[LIFE_TABLE];

void init_life( ) {
	srand(0);

	for( int i = 0; i < LIFE_TABLE; i++ ) {
		int t = __builtin_popcount(i);
		if( t == 3 )
			life_step[i] = true;
		if( t == 4 )
			life_step[i] = (i & (1 << 4));
	}
}

/*
  Life evolution, a pretty key subroutine
 */
template <int X, int Y>
grid<X,Y> evolve_once( const grid<X,Y>& start ) {
	grid<X,Y> stop;
	stop.cx = start.cx;
	stop.cy = start.cy;
	stop.bx = start.bx;
	stop.by = start.by;

	for( int x = 0; x < X; x++ ) {
	for( int y = 0; y < Y; y++ ) {
		int t = 0;
		for( int dx = -1; dx <= 1; dx++ ) {
		for( int dy = -1; dy <= 1; dy++ ) {
			t <<= 1;
			t += start.get_bool_uncentered( x+dx, y+dy );
		}
		}
		stop.set_uncentered( x, y, life_step[t] );
	}
	}

	return stop;
}
/*
  Evolve a position k times and record the results in an
  already-allocated array of size k+1.  Index 0 corresponds to the
  original position, which has already been assigned.

  (We need this for a specific case later, so as to avoid copying and
  memory allocations.)
*/
template <int X, int Y>
void evolve_many_special( grid<X,Y>* dest, const int& k ) {
	for( int i = 0; i < k; i++ ) {
		dest[i+1] = evolve_once(dest[i]);
	}
}

#endif
