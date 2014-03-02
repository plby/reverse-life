all: bin/main

bin:
	mkdir bin

bin/main: src/main.cc src/*.hh bin
	g++ -Wall -g -std=c++0x -isystem src/minisat/ src/main.cc src/minisat/core/Solver.or -o bin/main

bin/maino: src/main.cc src/*.hh bin/main
	g++ -O3 -std=c++0x -isystem src/minisat/ src/main.cc src/minisat/core/Solver.or -o bin/maino
