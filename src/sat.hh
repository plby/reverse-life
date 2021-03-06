#ifndef SAT_HH
#define SAT_HH

#include "core/Solver.h"
using namespace Minisat;
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
		cout << S.nClauses() << "\n";
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
	Lit FALSE;
	Solver S;

	// A list of known solutions:
	vector<vector<bool> > solutions;

	// This is a somewhat tricky construction for keeping track
	// whether things are possible
//	bool* possible;

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
		FALSE = get_literal( DEPTH+1, 0, 0 );
		S.addClause( ~FALSE );
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
					Lit temp = get_literal( depth+1, x+dx, y+dy );
					prev.push_back( temp );
				}
				}
				Lit next = get_literal( depth, x, y );
				include_clauses(S, prev, next);
			}
			}
		}
	}

	int get_numbering( int depth, int x, int y ) {
		return (depth * X + x) * Y + y;
	}
	int get_numbering_safe( int depth, int x, int y ) {
		if( 0 <= x and x < X and 0 <= y and y < Y ) {
			return (depth * X + x) * Y + y;
		} else {
			return -1;
		}
	}

	Lit get_literal( int depth, int x, int y, bool sign = false ) {
		if( 0 <= x and x < X and 0 <= y and y < Y ) {
			return mkLit( get_numbering( depth, x, y ), sign );
		} else {
			return FALSE;
		}
	}

	void solve( ) {
		// if( not S.simplify() ) {
		// 	cerr << "Did not expect unsatisfiable problem.\n";
		// 	exit( 50 );
		// }

		cout << "p cnf " << S.nVars() << " " << S.nClauses() << "\n";

		bool satisfiable = S.solve();
		if( not satisfiable ) {
			cerr << "Did not expect unsatisfiable problem.\n";
			exit( 60 );
		}

//		bool indeterminate = false;
		vector<bool> solution;
		for( int i = 0; i < S.nVars(); i++ ) {
			lbool l = S.modelValue( mkLit(i) );
			if( l == l_True ) {
				solution.push_back( true  );
			} else if( l == l_False ) {
				solution.push_back( false );
			} else {
				solution.push_back( false );
//				indeterminate = true;
			}
		}
		solutions.push_back(solution);
	}

	void extract_grid( grid<X,Y>& result ) {
		for( int x = 0; x < X; x++ ) {
		for( int y = 0; y < Y; y++ ) {
			result.set_safe_uncentered( x, y, solutions[0][get_numbering(DELTA,x,y)] );
		}
		}
	}
};

void sat( ) {
	const int WIDTH  = 2;
	const int HEIGHT = 2;
	const int MANY = 1;
	for( encoding e = 0; e < (1 << (WIDTH*HEIGHT)); e++ ) {
		grid<WIDTH,HEIGHT> g[MANY+1];
		g[0] = decode<WIDTH,HEIGHT>( e );
		for( int i = 0; i < MANY; i++ ) {
			g[i+1] = evolve_once<WIDTH,HEIGHT>( g[i] );
		}

		cout << g[MANY];
		life_solver<WIDTH,HEIGHT> ls( g[MANY], MANY, MANY );
		ls.solve();

		grid<WIDTH,HEIGHT> h[MANY+1];
		ls.extract_grid( h[0] );
		for( int i = 0; i < MANY; i++ ) {
			h[i+1] = evolve_once<WIDTH,HEIGHT>( h[i] );
		}

		if( g[MANY] != h[MANY] ) {
			cerr << "FAILED!!!\n";
			cerr << "e = " << e << "\n";
			for( int i = 0; i <= MANY; i++ ) {
				cerr << "g[" << i << "]=" << g[i];
			}
			for( int i = 0; i <= MANY; i++ ) {
				cerr << "h[" << i << "]=" << h[i];
			}
			exit( 20 );
		}
	}
}

#endif
