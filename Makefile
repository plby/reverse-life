all: bin/main bin/maino data/brain

bin:
	mkdir bin

bin/main: src/main.cc bin
	g++ -Wall -g -std=c++0x src/main.cc -o bin/main

bin/maino: src/main.cc bin/main
	g++ -O3 -std=c++0x src/main.cc -o bin/maino

data/brain:
	time head -c 1527353960 < /dev/zero > data/brain
