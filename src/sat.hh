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
	int DEPTH;   // how many steps back are we going
	int HORIZON; // how many steps back do we care about forcing
	Lit TRUE;
	Solver S;

	// A list of known solutions:
	vector<vector<bool> > solutions;

	// This is a somewhat tricky construction for keeping track
	// whether things are possible
	bool* possible;

	life_solver( grid<X,Y> final, int DEPTH, int HORIZON ) :
		final(final), DEPTH(DEPTH), HORIZON(HORIZON) {
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
		assert( S.nVars() == (DEPTH+1) * X * Y + 1 );

		// Start record-keeping about possible variables
//		possible = vector<int>( DEPTH * X * Y + 1, 0 );

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

	int get_numbering( int depth, int x, int y ) {
		return (depth * X + x) * Y + y;
	}

	Lit get_literal( int depth, int x, int y, bool sign = false ) {
		if( 0 <= x and x < X and 0 <= y and y < Y ) {
			return mkLit( get_numbering( depth, x, y ), sign );
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

		bool indeterminate = false;
		vector<bool> solution;
		for( int i = 0; i < S.nVars(); i++ ) {
			lbool l = S.modelValue( mkLit(i) );
			if( l == l_True ) {
				solution.push_back( true  );
			} else if( l == l_False ) {
				solution.push_back( false );
			} else {
				solution.push_back( false );
				indeterminate = true;
			}
		}
		solutions.push_back(solution);
	}

	void extract_grid( grid<X,Y>& result ) {
		for( int x = 0; x < X; x++ ) {
		for( int y = 0; y < Y; y++ ) {
			result.set_uncentered( x, y, solutions[0][get_numbering(1,x,y)] );
		}
		}
	}
};

void sat( ) {
	grid<4,4> g = decode<4,4>( uniform_smallint( 1 << 16 ) );
	life_solver<4,4> ls( g, 1, 1 );
	ls.solve();
	grid<4,4> h;
	ls.extract_grid( h );
	cout << g << "\n" << h << "\n";
}

#endif
