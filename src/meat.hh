#ifndef MEAT_HH
#define MEAT_HH

#include "grid.hh"

/*
This is where the real meat of the predictor is.
*/

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

	brain_data( string file ) {
		if( DISABLE_BRAIN )
			return;

		int fd;
		if( (fd = open(file.c_str(), O_RDWR | O_CREAT, S_IRUSR | S_IWUSR)) < 0 ) {
			cerr << "Couldn't open " << file << " file.\n";
			exit( 1 );
		}

		struct stat st;
		if( fstat(fd, &st) < 0 ) {
			cerr << "Couldn't stat " << file << " file.\n";
			exit( 2 );
		}

		int want_len = BRAIN * sizeof(unsigned int);
		int have_len = st.st_size;
		if( have_len != want_len ) {
			if( have_len > 0 ) {
				cerr << "File " << file << " is neither the expected size nor empty.\n";
				cerr << "Expected " << want_len << "  actual " << have_len << "  (bytes)\n";
				exit( 10 );
			}
			if( ftruncate(fd, want_len) != 0 ) {
				cerr << "Couldn't set " << file << " file size.\n";
				exit( 3 );
			}
		}

		data = (unsigned int *)mmap( 0, want_len,
					     PROT_READ | PROT_WRITE,
					     MAP_FILE  | MAP_SHARED, fd, 0 );
		if( data == MAP_FAILED ) {
			cerr << "Couldn't mmap " << file << " file.\n";
			exit( 4 );
		}
		if( close(fd) != 0 ) {
			cerr << "Couldn't close " << file << " file.\n";
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
		if( DISABLE_BRAIN )
			return;
		if( munmap(data, BRAIN * sizeof(unsigned int)) == -1 ) {
			cerr << "Error in munmap of a brain.\n";
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
} brain( "data/brain" );

void train_once( training_data<N,N> d ) {
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
		        big_grid h = d.gs[BURN+delta]; 
			grid<K,K> g = h.subgrid<K,K>( x, y, 2, 2 );
			encoding e = encode<K,K>( g );

			brain.add( delta, bucket(d.p), e, truth );
		}
	}
	}
}
void train_many( int TRAIN = 50000, int TRAIN_REPORT = 10000 ) {
	for( int i = 0; i < TRAIN; i++ ) {
		if( TRAIN_REPORT > 0 and (i % TRAIN_REPORT) == 0 )
			cout << i << " training grids generated." << endl;
		training_data<N,N> d;
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
bool predict_from_bucket( int delta, int bucket, grid<K,K> g ) {
	return predict_from_bucket( delta, bucket, encode<K,K>(g) );
}

const double TOLERANCE = 0.1;
bool predict_with_likelihood( int delta, big_grid stop, int x, int y,
			      double likelihood[BUCKETS], bool no_recurse );
double hack( int delta, big_grid stop, int x, int y, double likelihood[BUCKETS] ) {
	int t = 0;
	for( int dx = -1; dx <= 1; dx++ ) {
	for( int dy = -1; dy <= 1; dy++ ) {
		t *= 2;
		t += stop.represented_uncentered( x+dx, y+dy ) and
			predict_with_likelihood( delta+1, stop, x, y, likelihood, true );
	}
	}
	return life_step[t];
}

bool predict_with_likelihood( int delta, big_grid stop, int x, int y,
			      double likelihood[BUCKETS], bool no_recurse = false ) {
	grid<K,K> g = stop.subgrid<K,K>( x, y, 2, 2 );
	encoding e = encode<K,K>( g );

	double p = 0;
	for( int i = 0; i < BUCKETS; i++ ) {
		p += likelihood[i] * p_alive_from_bucket( delta, i, e );
	}

	if( (not no_recurse) and abs(p-0.5) < TOLERANCE and delta < 5 ) {
		return hack( delta, stop, x, y, likelihood );
	}

	return p > 0.5;
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
		result.set_uncentered(
			x, y,
			predict_with_likelihood( delta, stop, x, y, likelihood )
			);
	}
	}
	return result;
}

#endif
