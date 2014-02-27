all: binaries data/brain

binaries: 
	cd .. && $(MAKE)

data/brain:
	time head -c 1527353960 < /dev/zero > data/brain
