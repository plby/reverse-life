#ifndef ACTIONS_HH
#define ACTIONS_HH

void test( int TEST = 50000 ) {
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

void submit( predictor p, int SUBMIT = 50000 ) {
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

#endif
