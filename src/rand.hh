#ifndef RAND_HH
#define RAND_HH

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

#endif
