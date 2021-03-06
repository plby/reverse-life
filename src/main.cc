#include <cmath>
#include <bitset>
#include <vector>
#include <cassert>
#include <climits>
#include <cstdint>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <algorithm>
using namespace std;

/* mmap! */
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

// My headers
const bool DISABLE_BRAIN = false;

#include "grid.hh"
#include "life.hh"
#include "code.hh"
#include "rand.hh"
#include "generate.hh"
#include "consts.hh"
#include "predictors.hh"
#include "meat.hh"
#include "actions.hh"

void init( ) {
	init_life();
	init_code();
}

int main( ) {
	init();

	train_many();
	test();

//	submit( predict );

	return 0;
}
