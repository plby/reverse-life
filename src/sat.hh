#ifndef SAT_HH
#define SAT_HH

#include "all.cc"
#include "life.hh"

void include_clauses( Solver& S, const vector<Lit>& prev, const Lit& next ) {
	for( int i = 0; i < LIFE_TABLE; i++ ) {
		vec<Lit> clause;
		for( int j = 0; j < 9; j++ ) {
			if( i & (1 << j) ) {
				clause.push( ~prev[j] );
			} else {
				clause.push(  prev[j] );
			}
		}
		if( life_step[i] ) {
			clause.push(  next );
		} else {
			clause.push( ~next );
		}
		S.addClause( clause );
	}
}

/*
  A single life_solver instance wraps a variety of functions involving
  using a SAT solver to invert game of life positions.
*/
template<int X, int Y>
struct life_solver {
	grid<X, Y> final; // last life configuration
	int DEPTH; // how many steps back are we going
	Lit TRUE;
	Solver S;

	// A list of known solutions:
	vector<vector<bool> > solutions;

	bool* possible;

	life_solver( grid<X,Y> f ) : final(f), DEPTH(0) {
		// Add one variable for each cell at each of DEPTH+1
		// times
		for( int depth = 0; depth <= DEPTH; depth++ ) {
			for( int x = 0; x < X; x++ ) {
			for( int y = 0; y < Y; y++ ) {
				S.newVar();
			}
			}
		}
		// Add a variable that is forced true, for convenience
		S.newVar();
		TRUE = get_literal( DEPTH+1, 0, 0 );
		S.addClause( TRUE );
		assert( S.nVars() == DEPTH * X * Y + 1 );

		// Start record-keeping about possible variables
		possible = vector<int>( DEPTH * X * Y + 1, 0 );

		// Set up SAT problem boundary condition
		for( int x = 0; x < X; x++ ) {
		for( int y = 0; y < Y; y++ ) {
			if( final.get_bool_uncentered(x,y) ) {
				S.addClause( get_literal( 0, x, y ) );
			} else {
				S.addClause( get_literal( 0, x, y, true ) );
			}
		}
		}

		// Set up SAT problem transition table
		for( int depth = 0; depth < DEPTH; depth++ ) {
			for( int x = 0; x < X; x++ ) {
			for( int y = 0; y < Y; y++ ) {
				vector<Lit> prev;
				for( int dx = -1; dx <= +1; dx++ ) {
				for( int dy = -1; dy <= +1; dy++ ) {
					Lit temp = get_literal( depth, x+dx, y+dy );
					prev.push_back( temp );
				}
				}
				Lit next = get_literal( depth+1, x, y );
				include_clauses(S, prev, next);
			}
			}
		}
	}

	Lit get_literal( int depth, int x, int y, bool sign = false ) {
		if( 0 <= x and x < X and 0 <= y and y < Y ) {
			int t = (depth * X + x) * Y + y;
			return mkLit( t, sign );
		} else {
			return ~TRUE;
		}
	}

	vector<bool> solve( ) {
		vector<bool> result;

		if( !S.simplify() ) {
			cerr << "Did not expect unsatisfiable problem.\n";
			exit( 50 );
		}

		vec<Lit> dummy;
		bool satisfiable = S.solve();
		if( not satisfiable )
			return result;

		for( int i = 0; i < S.nVars(); i++ ) {
			lbool l = S.modelValue( mkLit(i) );
//			if( l == 
		}
		for( int i = 0; i < S.nVars(); i++) {
			lbool l = S.modelValue( mkLit(i) );
			printf( " " );
			if( l == l_True ) {
				printf( "+" );
			} else if( l == l_False ) {
				printf( "-" );
			} else {
				printf( "?" );
			}
			printf( "%d", i+1 );
		}
		
	}
};

void sat( ) {
}

#endif
