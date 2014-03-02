#ifndef GENERATE_H
#define GENERATE_H

#include "consts.h"

/*
  This data structure contains the information necessary for training,
  and a subset of it is useful for testing.
*/
const int GRIDS = BURN + DELTA + 1;

template <int X, int Y>
struct training_data {
	double p;
	grid<X,Y> gs[GRIDS];

	/*
	  Generate a random position according to the Bernoulli(p)
	  with p uniform model, and evolve it BURN+DELTA steps.

	  Empty grids are not filtered out at this point.
	 */
	training_data() {
		p = uniform_real( 0.01, 0.99 );
		for( int x = 0; x < X; x++ ) {
		for( int y = 0; y < Y; y++ ) {
			gs[0].set_uncentered( x, y, bernoulli(p) );
		}
		}
		evolve_many_special( gs, BURN+DELTA );
	}
};

template <int X, int Y>
struct testing_data {
	int delta;
	grid<X,Y> start;
	grid<X,Y> stop;

	/*
	  One reasonably important thing here is we filter out empty grids.
	*/
	testing_data( ) {
		do {
			training_data<X,Y> d;
			delta = uniform_smallint( 1, 5+1 ); // remember that the right end-point is excluded
			start = d.gs[BURN];
			stop  = d.gs[BURN+delta];
		} while( stop.g.count() == 0 );
	}
};

#endif
