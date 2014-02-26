all: main maino

main: main.cc
	g++ -Wall -g -std=c++11 main.cc -o main

maino: main.cc main
	g++ -O3 -std=c++11 main.cc -o maino
