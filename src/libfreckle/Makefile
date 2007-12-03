CPP=g++
CC=gcc
CPPFLAGS=-fPIC -Wall
LDFLAGS=-shared -Wl,-soname,test.so -lc

all: test.so

test.o: test.c
	$(CPP) $(CPPFLAGS) -c test.c

test.so: test.o
	$(CPP) $(CPPFLAGS) $(LDFLAGS) test.o -o test.so


DotStorageChunk.o: DotStorageChunk.cpp DotStorageChunk.h
	$(CPP) $(CPPFLAGS) -c DotStorageChunk.cpp


#
# CPP Test Suites
#
testDotStorageChunk.cpp: testDotStorageChunk.h DotStorageChunk.cpp DotStorageChunk.h
	./cxxtestgen.pl --error-printer -o testDotStorageChunk.cpp testDotStorageChunk.h

testDotStorageChunk: testDotStorageChunk.cpp DotStorageChunk.o
	$(CPP) -I./ -o testDotStorageChunk testDotStorageChunk.cpp DotStorageChunk.o

runtests: testDotStorageChunk
	./testDotStorageChunk



#
# Cleans
#
cleantests:
	-rm testDotStorageChunk.cpp testDotStorageChunk


clean:
	-rm *.so
	-rm *.o

#gcc -fPIC -g -c -Wall a.c
#gcc -fPIC -g -c -Wall b.c
#gcc -shared -Wl,-soname,libmystuff.so.1 \
#    -o libmystuff.so.1.0.1 a.o b.o -lc