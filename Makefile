all: bin/main bin/maino

bin/main: src/main.cc
	g++ -Wall -g -std=c++11 src/main.cc -o bin/main

bin/maino: src/main.cc bin/main
	g++ -O3 -std=c++11 src/main.cc -o bin/maino
