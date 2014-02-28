all: bin/main bin/maino data/brain

bin:
	mkdir bin

bin/main: src/main.cc bin
	g++ -Wall -g -std=c++0x src/main.cc -o bin/main

bin/maino: src/main.cc bin/main
	g++ -O3 -std=c++0x src/main.cc -o bin/maino

data/brain:
	mkdir -p data/ && time head -c 6109415840 < /dev/zero > data/brain
