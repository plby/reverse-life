all: bin/main

bin:
	mkdir bin

bin/main: src/main.cc src/*.hh bin
	g++ -Wall -g -std=c++0x -isystem src/glucose/ src/main.cc -o bin/main

bin/maino: src/main.cc src/*.hh bin/main
	g++ -O3 -std=c++0x -isystem src/glucose/ src/main.cc -o bin/maino
