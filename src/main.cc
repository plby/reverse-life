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

// My headers
#include "grid.h"
#include "life.h"
#include "code.h"
#include "rand.h"
#include "generate.h"
#include "consts.h"

const int SUBMIT = 50000;
const int TEST   = 50000;
const int TRAIN  = 50000;

const int TEST_REPORT  = 10000;
const int TRAIN_REPORT = 10000;

/*
  C++ doesn't have first-class function objects, and I don't want to
  wrap this into a class, so a predictor is just a typedef to a
  function pointer.
*/
typedef grid<N,N> big_grid;
typedef big_grid (*predictor)( int, big_grid );

/*
  Return the number of incorrect cells guessed
*/
template <int X, int Y>
int grade_once( testing_data<X,Y> d, big_grid guess ) {
	int result = 0;
	for( int x = 0; x < X; x++ ) {
	for( int y = 0; y < Y; y++ ) {
		result += (d.start.get_bool_uncentered(x,y)
			   !=guess.get_bool_uncentered(x,y) );
	}
	}
	return result;
}
template <int X, int Y>
int grade_once( testing_data<X,Y> d, predictor p ) {
	big_grid guess = p( d.delta, d.stop );
	return grade_once<X,Y>( d, guess );
}

/*
  Grade some prediction functions
*/
template <int X, int Y>
vector<double> grade_many( vector<predictor> ps, int trials = 100000 ) {
	int P = ps.size();
	vector<int> wrong( P );
	vector<int> total( P );
	for( int i = 0; i < trials; i++ ) {
		testing_data<X,Y> d;
		for( int j = 0; j < P; j++ ) {
			wrong[j] += grade_once( d, ps[j] );
			total[j] += X*Y;
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
	unsigned int& get( const int& delta, const int& bucket, const grid<K,K>& grid, const bool entry ) {
		encoding code = encode<K,K>( grid );
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

template <int X, int Y>
void train_once( training_data<X,Y> d ) {
	// Find maximum delta without dead grid
	int alive;
	for( alive = 1; alive <= DELTA; alive++ ) {
		if( d.gs[BURN+alive].g.count() == 0 )
			break;
	}

	for( int x = 0; x < X; x++ ) {
	for( int y = 0; y < Y; y++ ) {
		bool truth = d.gs[BURN].get_bool_uncentered( x, y );

		for( int delta = 1; delta < alive; delta++ ) {
			grid<K,K> g = d.gs[BURN+delta].subgrid<K,K>( x, y, 2, 2 );
			encoding e = encode<K,K>( g );

			brain.add( delta, bucket(d.p), e, truth );
		}
	}
	}
}
template <int X, int Y>
void train_many( ) {
	for( int i = 0; i < TRAIN; i++ ) {
		if( TRAIN_REPORT > 0 and (i % TRAIN_REPORT) == 0 )
			cout << i << " training grids generated." << endl;
		training_data<X,Y> d;
		train_once<X,Y>( d );
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
bool predict_from_bucket( int delta, int bucket, grid<K,K> g ) {
	return predict_from_bucket( delta, bucket, encode<K,K>(g) );
}

big_grid predict( int delta, big_grid stop ) {
	/*
	  Predict the bucket for p first using naive Bayes.
	*/
	double log_likelihood[BUCKETS];
	for( int i = 0; i < BUCKETS; i++ )
		log_likelihood[i] = 0;

	for( int x = 0; x < N; x++ ) {
	for( int y = 0; y < N; y++ ) {
		grid<K,K> g = stop.subgrid<K,K>( x, y, 2, 2 );
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
		grid<K,K> g = stop.subgrid<K,K>( x, y, 2, 2 );
		encoding e = encode<K,K>( g );

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

	vector<double> result = grade_many<N,N>( ps, TEST );
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
		for( int i = 1; i <= N*N; i++ )
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
		if( not( 1 <= delta and delta <= DELTA ) )
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

//	test();

//	submit( predict );

	return 0;
}
