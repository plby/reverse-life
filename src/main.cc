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
#include "grid.h"
#include "life.h"
#include "code.h"
#include "rand.h"
#include "generate.h"
#include "consts.h"
#include "predictors.h"
#include "meat.h"
#include "actions.h"

void init( ) {
	init_life();
	init_code();
}

int main( ) {
	init();

//	test();

//	submit( predict );

	return 0;
}
