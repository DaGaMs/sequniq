CPP = clang++
CFLAGS = -arch x86_64 -std=gnu++11 -stdlib=libc++ -Os -fasm-blocks -fstrict-aliasing -fvisibility=hidden -fvisibility-inlines-hidden
INCLUDE = -IExternal/TCLAP/include -IExternal/kseq

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