#ifndef PREDICTORS_HH
#define PREDICTORS_HH

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
vector<double> grade_many( vector<predictor> ps, int TEST = 50000, int TEST_REPORT = 10000 ) {
	int P = ps.size();
	vector<int> wrong( P );
	vector<int> total( P );
	for( int i = 0; i < TEST; i++ ) {
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

#endif
