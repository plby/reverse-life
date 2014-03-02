#ifndef CODE_HH
#define CODE_HH

#include "grid.hh"

/*
  The main encoding is fairly simple.  It does not encode (cx,cy).  So
  it just lists, in order, the different possible (bx,by,grid) tuples,
  beginning with (bx,by) = (0,0) and afterwards proceeding in
  lexicographic order -- eg, for a centered 5x5 grid, beginning with
  (-2,-2) and (-2,-1).  See the code for more details.
 */
typedef uint32_t encoding;
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

#endif
