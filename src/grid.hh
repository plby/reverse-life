#ifndef GRID_HH
#define GRID_HH

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

#endif
