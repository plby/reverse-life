#include <cmath>
#include <bitset>
#include <vector>
#include <cassert>
#include <climits>
#include <cstdint>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <algorithm>
using namespace std;

/* mmap! */
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

const int N = 20;
const int M = N * N;

const int SUBMIT = 50000;
const int TEST   = 50000;
const int TRAIN  = 50000;

const int TEST_REPORT  = 10000;
const int TRAIN_REPORT = 10000;

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

	void set_safe_uncentered( const int& x, const int& y, const bool& v ) {
		if( inside_uncentered(x,y) )
			g.set( x*Y + y, v );
	}
	void set_uncentered( const int& x, const int& y, const bool& v ) {
		g.set( x*Y + y, v );
	}
	void set( const int& x, const int& y, const bool& v ) {
		set_uncentered( x+cx, y+cy, v );
	}

	/*
	  Returns a smaller grid, centered at (x,y) with the
	  underlying grid offset by (cx2, cy2).
	*/
	template<int X2, int Y2>
	grid<X2, Y2> subgrid( const int& x, const int& y, const int& cx2, const int& cy2 ) {
		grid<X2, Y2> result;
		result.cx = cx2;
		result.cy = cy2;

		// We make some assumptions
		assert( bx == 0 and by == 0 );

		// The computation of the subordinate (bx, by) is somewhat tricky
		int t;
		{
			t = (x - cx2);
			if( t < 0 )
				result.bx = -t;
			t = (x - cx2) + (X2 - X);
			if( t > 0 )
				result.bx = -t;
		}
		{
			t = (y - cy2);
			if( t < 0 )
				result.by = -t;
			t = (y - cy2) + (Y2-Y);
			if( t > 0 )
				result.by = -t;
		}
		if( abs(result.bx) >= X or abs(result.by) >= Y ) {
			result.bx = X;
			result.by = Y;
		}

		// Copy over data
		for( int dx = 0; dx < X; dx++ ) {
		for( int dy = 0; dy < Y; dy++ ) {
			result.set_safe_uncentered( dx, dy,
						    get_bool_uncentered( x-cx2+dx,
									 y-cy2+dy )
				);
		}
		}

		return result;
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
bool operator == ( const grid<X,Y>& g, const grid<X,Y>& h ) {
	return not( g != h);
}
template <int X, int Y>
bool operator < ( const grid<X,Y>& g, const grid<X,Y>& h ) {
	if( g.bx != h.bx )
		return g.bx < h.bx;
	if( g.by != h.by )
		return g.by < h.by;
	
	for( int x = 0; x < X; x++ ) {
	for( int y = 0; y < Y; y++ ) {
		if( g.get_int_uncentered(x,y) != h.get_int_uncentered(x,y) ) {
			return g.get_int_uncentered(x,y) < h.get_int_uncentered(x,y);
		}
	}
	}
	return false;
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
  Very simple function to remove duplicates from a vector.
*/
template<typename t>
void sort_unique( vector<t>& v ) {
	sort( v.begin(), v.end() );
	v.erase( unique( v.begin(), v.end() ), v.end() );
}

/*
  Return all of the symmetrical versions of a given grid.
*/
template <int X, int Y>
grid<Y,X> flip_grid( const grid<X,Y>& in ) {
	grid<Y,X> out;
	out.by = in.bx;
	out.bx = in.by;
	out.cy = in.cx;
	out.cx = in.cy;
	if( in.bx == X and in.by == Y ) {
		return out;
	}
	for( int x = 0; x < X; x++ ) {
	for( int y = 0; y < Y; y++ ) {
		out.set_uncentered( y, x, in.get_bool_uncentered( x, y ) );
	}
	}
	return out;
}

template <int X, int Y>
grid<Y,X> rotate_grid( const grid<X,Y>& in ) {
	grid<Y,X> out;
	out.by = +in.bx;
	out.bx = -in.by;
	out.cy = in.cx;
	out.cx = (Y-1)-in.cy;
	if( in.bx == X and in.by == Y ) {
		out.bx = in.by;
		return out;
	}
	for( int x = 0; x < X; x++ ) {
	for( int y = 0; y < Y; y++ ) {
		out.set_uncentered( (Y-1)-y, x, in.get_bool_uncentered( x, y ) );
	}
	}
	return out;
}

template <int X>
vector<grid<X,X> > symmetric( const grid<X,X>& first ) {
	vector<grid<X,X> > result;

	grid<X,X> temp = first;
	for( int i = 0; i < 2; i++ ) {
		result.push_back( temp );
		for( int j = 0; j < 3; j++ ) { // the 3 is a super small optimization
			temp = rotate_grid( temp );
			result.push_back( temp );
		}
		if( i == 0 )
			temp = flip_grid( temp );
	}

	sort_unique( result );
	return result;
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
  Take advantage of symmetries by only encoding positions up to
  symmetry.
*/
const int BIGCODE   = 38183849;
const int SMALLCODE =  4793312;
encoding *smallcode;
void init_code( ) {
	int fd;
	if( (fd = open("data/codemap", O_RDWR | O_CREAT, S_IRUSR | S_IWUSR)) < 0 ) {
		cerr << "Couldn't open data/codemap file.\n";
		exit( 1 );
	}

	struct stat st;
	if( fstat(fd, &st) < 0 ) {
		cerr << "Couldn't stat data/codemap file.\n";
		exit( 2 );
	}

	int want_len = BIGCODE*sizeof(encoding);
	int have_len = st.st_size;
	if( have_len != want_len ) {
		if( ftruncate(fd, want_len) != 0 ) {
			cerr << "Couldn't set data/codemap file size.\n";
			exit( 3 );
		}
	}

	smallcode = (encoding *)mmap( 0, want_len,
				      PROT_READ | PROT_WRITE,
				      MAP_FILE  | MAP_SHARED, fd, 0 );
	if( smallcode == MAP_FAILED ) {
		cerr << "Couldn't mmap data/codemap file.\n";
		exit( 4 );
	}
	if( close(fd) != 0 ) {
		cerr << "Couldn't close data/codemap file.\n";
		exit( 5 );
	}

	if( smallcode[BIGCODE-1] == SMALLCODE-1 ) {
		// If this is correct, then we will assume the rest of
		// the mapping is correct to save initialization time
		return;
	}

	// Otherwise redo the mapping
	encoding e, k;
	k = 1;
	for( e = 0; e < (encoding)BIGCODE; e++ ) {
		grid<5,5> g = decode<5,5>( e );
		vector<grid<5,5> > syms = symmetric( g );

		if( smallcode[encode<5,5>( syms[0] )] != 0 )
			continue;

		for( int i = 0; i < (int)syms.size(); i++ ) {
			encoding f = encode<5,5>( syms[i] );
			assert( smallcode[f] == 0 );
			smallcode[f] = k;
		}
		k++;
	}
	assert( k == (int)SMALLCODE );
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
	  One reasonably important thing here is we filter out empty grids.
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
		testing_data d;
		for( int j = 0; j < P; j++ ) {
			wrong[j] += grade_once( d, ps[j] );
			total[j] += M;
		}

		if( TEST_REPORT > 0 and i > 0 and (i % TEST_REPORT) == 0 ) {
			for( int j = 0; j < P; j++ ) {
				cout << setprecision(5)
				     << (double)(wrong[j]) / (double)(total[j])
				     << "\t";
			}
			cout << "(" << i << ")\n";
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

/*
  This is the first non-trivial predictor.  It uses mmap to manage its
  data representation.
*/
const int BUCKETS = 8;
const int ENTRIES = 2;
const int BRAIN = DELTA * BUCKETS * SMALLCODE * ENTRIES;

int bucket( double p ) {
	double cutoffs[BUCKETS] = 
	{
		0.164699,
		0.258496,
		0.352168,
		0.445686,
		0.539146,
		0.632559,
		0.728846,
		1
	};
	int result;
	for( result = 0;
	     result < BUCKETS-1 and p > cutoffs[result];
	     result++ ) {
		// empty loop
	}
	return result;
}

struct brain_data {
	unsigned int *data;

	brain_data( ) {
		int fd;
		if( (fd = open("data/brain", O_RDWR | O_CREAT, S_IRUSR | S_IWUSR)) < 0 ) {
			cerr << "Couldn't open data/brain file.\n";
			exit( 1 );
		}

		struct stat st;
		if( fstat(fd, &st) < 0 ) {
			cerr << "Couldn't stat data/brain file.\n";
			exit( 2 );
		}

		int want_len = BRAIN * sizeof(unsigned int);
		int have_len = st.st_size;
		if( have_len != want_len ) {
			if( have_len > 0 ) {
				cerr << "File data/brain is neither the expected size nor empty.\n";
				cerr << "Expected " << want_len << "  actual " << have_len << "  (bytes)\n";
				exit( 10 );
			}
			if( ftruncate(fd, want_len) != 0 ) {
				cerr << "Couldn't set data/brain file size.\n";
				exit( 3 );
			}
		}

		data = (unsigned int *)mmap( 0, want_len,
					     PROT_READ | PROT_WRITE,
					     MAP_FILE  | MAP_SHARED, fd, 0 );
		if( data == MAP_FAILED ) {
			cerr << "Couldn't mmap data/brain file.\n";
			exit( 4 );
		}
		if( close(fd) != 0 ) {
			cerr << "Couldn't close data/brain file.\n";
			exit( 5 );
		}

		// dummy read
		unsigned long long int entries = 0;
		for( int i = 0; i < BRAIN; i++ ) {
			entries += data[i];
		}
		cerr << entries << " entries currently in the brain.\n";
	}

	~brain_data( ) {
		if( munmap(data, BRAIN * sizeof(unsigned int)) == -1 ) {
			cerr << "Error in munmap of data/brain.\n";
			exit( 7 );
		}
	}

	unsigned int& get( const int& delta, const int& bucket, const encoding& code, const bool entry ) {
		int index = (((delta-1) * BUCKETS + bucket) * SMALLCODE + smallcode[code]) * ENTRIES + entry;
		return data[index];
	}
	unsigned int& get( const int& delta, const int& bucket, const grid<5,5>& grid, const bool entry ) {
		encoding code = encode<5,5>( grid );
		return get( delta, bucket, code, entry );
	}
	void add( const int& delta, const int& bucket, const encoding& code, const bool entry ) {
		unsigned int& g = get( delta, bucket, code, entry );
		if( g > UINT_MAX - 10 ) {
			cerr << "Brain entry is too big.\n";
			exit( 1 );
		}
		g++;
	}
} brain;

void train_once( training_data d ) {
	// Find maximum delta without dead grid
	int alive;
	for( alive = 1; alive <= DELTA; alive++ ) {
		if( d.gs[BURN+alive].g.count() == 0 )
			break;
	}

	for( int x = 0; x < N; x++ ) {
	for( int y = 0; y < N; y++ ) {
		bool truth = d.gs[BURN].get_bool_uncentered( x, y );

		for( int delta = 1; delta < alive; delta++ ) {
			grid<5,5> g = d.gs[BURN+delta].subgrid<5,5>( x, y, 2, 2 );
			encoding e = encode<5,5>( g );

			brain.add( delta, bucket(d.p), e, truth );
		}
	}
	}
}
void train_many( ) {
	for( int i = 0; i < TRAIN; i++ ) {
		if( TRAIN_REPORT > 0 and (i % TRAIN_REPORT) == 0 )
			cout << i << " training grids generated." << endl;
		training_data d;
		train_once( d );
	}
}

double p_alive_from_bucket( int delta, int bucket, encoding e ) {
	int dead, alive;
	if( 0 <= bucket and bucket < BUCKETS ) {
		dead  = brain.get( delta, bucket, e, false );
		alive = brain.get( delta, bucket, e, true  );
	} else {
		dead  = 0;
		alive = 0;
		for( int i = 0; i < BUCKETS; i++ ) {
			dead  += brain.get( delta, i, e, false );
			alive += brain.get( delta, i, e, true  );			
		}
	}

	// The following reflects a minimal 1/7 prior probability of being dead
	dead += 5;

	return double(alive + 1) / double(dead + alive + 2);
}
bool predict_from_bucket( int delta, int bucket, encoding e ) {
	return p_alive_from_bucket( delta, bucket, e ) > 0.5;
}
bool predict_from_bucket( int delta, int bucket, grid<5,5> g ) {
	return predict_from_bucket( delta, bucket, encode<5,5>(g) );
}

/*
  The adjustment factor is not chosen too scientifically. :(
*/
const double ADJUST = 25;
big_grid predict( int delta, big_grid stop ) {
	/*
	  Predict the bucket for p first using naive Bayes.
	*/
	double log_likelihood[BUCKETS];
	for( int i = 0; i < BUCKETS; i++ )
		log_likelihood[i] = 0;

	for( int x = 0; x < N; x++ ) {
	for( int y = 0; y < N; y++ ) {
		grid<5,5> g = stop.subgrid<5,5>( x, y, 2, 2 );
		for( int i = 0; i < BUCKETS; i++ ) {
			int dead  = brain.get( delta, i, g, false );
			int alive = brain.get( delta, i, g, true  );
			int total = dead + alive;
			log_likelihood[i] += log(1 + total) / ADJUST;
		}
	}
	}	

	// Normalize a bit to keep things in a manageable range
	double biggest = 0;
	for( int i = 0; i < BUCKETS; i++ ) {
		if( biggest < log_likelihood[i] )
			biggest = log_likelihood[i];
	}
	for( int i = 0; i < BUCKETS; i++ ) {
		log_likelihood[i] -= biggest;
	}
	double likelihood[BUCKETS];
	for( int i = 0; i < BUCKETS; i++ ) {
		likelihood[i] = exp(log_likelihood[i]);
	}
	double sum = 0;
	for( int i = 0; i < BUCKETS; i++ ) {
		sum += likelihood[i];
	}
	for( int i = 0; i < BUCKETS; i++ ) {
		likelihood[i] /= sum;
	}

	/*
	  Now use that information to make a better prediction for
	  each cell.
	*/
	big_grid result;
	for( int x = 0; x < N; x++ ) {
	for( int y = 0; y < N; y++ ) {
		grid<5,5> g = stop.subgrid<5,5>( x, y, 2, 2 );
		encoding e = encode<5,5>( g );

		double p = 0;
		for( int i = 0; i < BUCKETS; i++ ) {
			p += likelihood[i] * p_alive_from_bucket( delta, i, e );
		}

		result.set_uncentered( x, y, p > 0.5 );
	}
	}
	return result;
}

void test( ) {
	vector<predictor> ps;

	cout << "stop\t";
	ps.push_back( start_at_stop );

	cout << "dead\t";
	ps.push_back( all_dead );

	cout << "predict\t";
	ps.push_back( predict );

	cout << "\n";

	vector<double> result = grade_many( ps, TEST );
	for( int i = 0; i < (int)result.size(); i++ ) {
		cout << setprecision(5) << result[i] << "\t";
	}
	cout << "\n";
}

void fail_parse( string detail ) {
	cerr << "Failed to parse test file because of " << detail << ".\n";
	exit( 2 );
}

void getcomma( istream& in ) {
	if( in.get() != ',' ) {
		fail_parse( "comma" );
	}
}

void submit( predictor p ) {
	ifstream fin( "data/raw/test.csv" );
	// Ignore header line:
	string header;
	getline( fin, header );

	ofstream fout( "data/raw/submit.csv" );
	// Header line:
	{
		fout << "id";
		for( int i = 1; i <= M; i++ )
			fout << ",start." << i;
		fout << "\n";
	}

	for( int i = 1; i <= SUBMIT; i++ ) {
		// Read a single line
		int id;
		int delta;
		big_grid stop;

		fin >> id;
		if( i != id )
			fail_parse( "id" );
		getcomma(fin);
		fin >> delta;
		if( not( 1 <= delta and delta <= 5 ) )
			fail_parse( "delta" );
		for( int x = 0; x < N; x++ ) {
       		for( int y = 0; y < N; y++ ) {
			getcomma(fin);
			int t;
			fin >> t;
			if( t != 0 and t != 1 )
				fail_parse( "bool" );
			stop.set_uncentered( x, y, t );
		}
		}

		// Make a prediction
		big_grid guess = p( delta, stop );

		// Output
		fout << id;
		for( int x = 0; x < N; x++ ) {
       		for( int y = 0; y < N; y++ ) {
			fout << "," << guess.get_int_uncentered( x, y );
		}
		}
		fout << "\n";
	}
}

void init( ) {
	init_life();
	init_code();
}

int main( ) {
	init();

	test();

//	submit( predict );

	return 0;
}
