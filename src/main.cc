#include <bitset>
#include <vector>
#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>
using namespace std;

const int N = 20;
const int M = N * N;

const int CODE5 = 38183849;

/*
  We precompute a lookup table for the game of life function for
  speed's sake.
 */
const int LIFE_TABLE = 1 << 9;
bool life_step[LIFE_TABLE];

void init( ) {
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
  is represented.  Specifically, a positive bx indicates how many
  columns on the left are outside of the grid, while a negative the
  opposite on the other side.  A bx of 0 means at least part of each
  column is inside.  The same for by, where positive by means rows
  missing at the top.  Example:

  X....  X = outside
  X....  . = inside
  X....
  XXXXX  bx = +1  cx = 0
  XXXXX  by = -2  cy = 0

  One relatively important thing is that every 5x5 such grid can be
  transformed into a 26-bit representation, which will be used to
  index into an array.
 */

template <int X, int Y>
struct grid {
	/* (0,0) in our x-y coordinate system corresponds to (cx,cy)
	   in the usual system with (0,0) at the corner */
	int cx, cy;

	/* This represents how much of the grid is "out of bounds", in
	   the manner described above */
	int bx, by;

	bitset<X*Y> g;

	grid( ) : cx(0), cy(0), bx(0), by(0), g(0) {
	}

	bool represented_uncentered( const int& x, const int& y ) const {
		return 0 <= x and x < X and 0 <= y and y < Y;
	}
	bool represented( const int& x, const int& y ) const {
		return represented_uncentered( x+cx, y+cy );
	}

	bool inside_uncentered( const int& x, const int& y ) const {
		return represented_uncentered(x,y) and
			0 <= (x-bx) and (x-bx) < X and 0 <= (y-by) and (y-by) < Y;
	}
	bool inside( const int& x, const int& y ) const {
		return inside_uncentered(x+cx, y+cy);
	}

	int get_int_uncentered( const int& x, const int& y ) const {
		if( not inside_uncentered(x,y) )
			return -1;
		return get_raw_uncentered(x,y);
	}
	int get_int( const int& x, const int& y ) const {
		return get_int_uncentered(x+cx, y+cy);
	}
	bool get_bool_uncentered( const int& x, const int& y ) const {
		return inside_uncentered(x,y) and get_raw_uncentered(x,y);
	}
	bool get_bool( const int& x, const int& y ) const {
		return get_bool_uncentered(x+cx, y+cy);
	}
	bool get_raw_uncentered( const int& x, const int& y ) const {
		return g[ x*Y + y ];
	}
	bool get_raw( const int& x, const int& y ) const {
		return get_raw_uncentered( x+cx, y+cy );
	}

	void set_uncentered( const int& x, const int& y, const bool& v ) {
		g.set( x*Y + y, v );
	}
	void set( const int& x, const int& y, const bool& v ) {
		set_uncentered( x+cx, y+cy, v );
	}
};
typedef grid<N,N> big_grid;
typedef uint32_t encoding;

template <int X, int Y>
bool operator != ( const grid<X,Y>& g, const grid<X,Y>& h ) {
	if( g.bx != h.bx )
		return true;
	if( g.by != h.by )
		return true;
	return g.g != h.g;
}

template <int X, int Y>
ostream& operator << ( ostream& out, const grid<X,Y>& g ) {
	// Header information
	out << "grid " << X << "x" << Y << " "
	    << "bxy=(" << g.bx << "," << g.by << ") "
	    << "cxy=(" << g.cx << "," << g.cy << ")\n";
	for( int y = 0; y < Y; y++ ) {
		for( int x = 0; x < X; x++ ) {
			int t = g.get_int_uncentered(x,y);
			out << "_.#"[t+1];
		}
		out << "\n";
	}

	return out;
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

/*
  The main encoding is fairly simple.  It does not encode (cx,cy).  So
  it just lists, in order, the different possible (bx,by,grid) tuples,
  beginning with (bx,by) = (0,0) and afterwards proceeding in
  lexicographic order -- eg, for a centered 5x5 grid, beginning with
  (-2,-2) and (-2,-1).  See the code for more details.
 */
template <int X, int Y>
encoding encode( const grid<X,Y>& g ) {
	// Optimized for bx = by = 0
	if( g.bx == 0 and g.by == 0 ) {
		return g.g.to_ulong();
	}

	// Otherwise, build up to number
	encoding e = (1 << (X*Y));
	for( int bx = -(X-1); bx <= (X-1); bx++ ) {
	for( int by = -(Y-1); by <= (Y-1); by++ ) {
		if( bx == 0 and by == 0 )
			continue;
		
		int width  = X - abs(bx);
		int height = Y - abs(by);
		int area   = width * height;
		if( bx != g.bx or by != g.by ) {
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
	assert( g.bx == X and g.by == Y );
	return e;
}
template <int X, int Y>
grid<X,Y> decode( const encoding& f ) {
	grid<X,Y> g;
	// Optimized for bx = by = 0
	if( f < (1 << (X*Y)) ) {
		g.bx = 0;
		g.by = 0;
		g.g = bitset<X*Y>(f);
		return g;
	}
	
	// Otherwise, build up to number
	encoding e = (1 << (X*Y));
	for( int bx = -(X-1); bx <= (X-1); bx++ ) {
	for( int by = -(Y-1); by <= (Y-1); by++ ) {
		if( bx == 0 and by == 0 )
			continue;
		
		int width  = X - abs(bx);
		int height = Y - abs(by);
		int area   = width * height;
		if( f >= e + (1 << area) ) {
			e += (1 << area);
		} else {
			g.bx = bx;
			g.by = by;

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
	// If we're here, then this better represent the empty grid
	assert( e == f );
	g.bx = X;
	g.by = Y;
	return g;
}

/*
  Routines to help with randomness.
 */
double uniform_real( ) {
	return (double)rand()/(double)RAND_MAX;
}
double uniform_real( double max ) {
	return max * uniform_real();
}
double uniform_real( double min, double max ) {
	return min + uniform_real( max-min );
}
/* max is always excluded */
int uniform_smallint( int max ) {
	return rand() % max; // may be a bad idea for large values of max
}
/* max is always excluded */
int uniform_smallint( int min, int max ) {
	return min + uniform_smallint(max-min);
}
bool bernoulli( double p ) {
	double t = uniform_real();
	return t < p;
}

/*
  This data structure contains the information necessary for training,
  and a subset of it is useful for testing.
*/
const int BURN  = 5; // number of steps to burn in each grid
const int DELTA = 5; // maximum value of delta
const int GRIDS = BURN + DELTA + 1;
struct training_data {
	double p;
	big_grid gs[GRIDS];

	/*
	  Generate a random position according to the Bernoulli(p)
	  with p uniform model, and evolve it BURN+DELTA steps.

	  Empty grids are not filtered out at this point.
	 */
	training_data() {
		p = uniform_real( 0.01, 0.99 );
		for( int x = 0; x < N; x++ ) {
		for( int y = 0; y < N; y++ ) {
			gs[0].set_uncentered( x, y, bernoulli(p) );
		}
		}
		evolve_many_special( gs, BURN+DELTA );
	}
};

struct testing_data {
	int delta;
	big_grid start;
	big_grid stop;

	/*
	  It is in this step and this step only that we filter out empty grids.
	*/
	testing_data( ) {
		do {
			training_data d;
			delta = uniform_smallint( 1, 5+1 ); // remember that the right end-point is excluded
			start = d.gs[BURN];
			stop  = d.gs[BURN+delta];
		} while( stop.g.count() == 0 );
	}
};

/*
  C++ doesn't have first-class function objects, and I don't want to
  wrap this into a class, so a predictor is just a typedef to a
  function pointer.
*/
typedef big_grid (*predictor)( int, big_grid );

/*
  Return the number of incorrect cells guessed
*/
int grade_once( testing_data d, big_grid guess ) {
	int result = 0;
	for( int x = 0; x < N; x++ ) {
	for( int y = 0; y < N; y++ ) {
		result += (d.start.get_bool_uncentered(x,y)
			   !=guess.get_bool_uncentered(x,y) );
	}
	}
	return result;
}
int grade_once( testing_data d, predictor p ) {
	big_grid guess = p( d.delta, d.stop );
	return grade_once( d, guess );
}

/*
  Grade some prediction functions
*/
vector<double> grade_many( vector<predictor> ps, int trials = 100000 ) {
	int P = ps.size();
	vector<int> wrong( P );
	vector<int> total( P );
	for( int i = 0; i < trials; i++ ) {
		for( int j = 0; j < P; j++ ) {
			testing_data d;
			wrong[j] += grade_once( d, ps[j] );
			total[j] += M;
		}
	}
	vector<double> result( P );
	for( int j = 0; j < P; j++ ) {
		result[j] = (double)(wrong[j]) / (double)(total[j]);
	}
	return result;
}

/*
  Benchmark predictors
*/
big_grid all_dead( int delta, big_grid stop ) {
	big_grid dead;
	return dead;
}

big_grid start_at_stop( int delta, big_grid stop ) {
	return stop;
}

big_grid all_alive( int delta, big_grid stop ) {
	big_grid alive;
	for( int x = 0; x < N; x++ ) {
	for( int y = 0; y < N; y++ ) {
		alive.set_uncentered(x,y,true);
	}
	}
	return alive;
}

int main( ) {
	init();

	vector<predictor> ps;

	cout << "dead\t";
	ps.push_back( all_dead );

	cout << "stop\t";
	ps.push_back( start_at_stop );

	cout << "alive\t";
	ps.push_back( all_alive );

	cout << "\n";
	while( 1 ) {
		vector<double> result = grade_many( ps, 100000 );
		for( int i = 0; i < (int)result.size(); i++ ) {
			cout << setprecision(5) << result[i] << "\t";
		}
		cout << endl;
	}

	return 0;
}
