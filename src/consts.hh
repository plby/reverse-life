#ifndef CONSTS_HH
#define CONSTS_HH

const int N = 20;
const int K = 5;

const int BURN  = 5; // number of steps to burn in each grid
const int DELTA = 5; // maximum value of delta

/*
  The adjustment factor is not chosen too scientifically. :(
*/
const double ADJUST = K*K;

#endif
