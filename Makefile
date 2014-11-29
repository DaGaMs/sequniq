CPP = c++
CFLAGS = -std=gnu++11 -Ofast
INCLUDE = -IExternal/TCLAP/include -IExternal/kseq

prefix=/usr/local

all: sequniq

sequniq.o: sequniq/MurmurHash3.h sequniq/sequniq.cpp
	mkdir -p build
	$(CPP) $(CFLAGS) $(INCLUDE) -c sequniq/sequniq.cpp -o build/sequniq.o

MurmurHash3.o: sequniq/MurmurHash3.h sequniq/MurmurHash3.cpp
	mkdir -p build
	$(CPP) $(CFLAGS) $(INCLUDE) -c sequniq/MurmurHash3.cpp -o build/MurmurHash3.o

sequniq: MurmurHash3.o sequniq.o
	mkdir -p bin
	$(CPP) -lz -o bin/sequniq build/sequniq.o build/MurmurHash3.o

clean:
	$(RM) -rf build bin/*
	
install: sequniq
	install -m 0755 sequniq $(prefix)/bin

.PHONY: install