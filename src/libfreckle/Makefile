CPP=g++
CC=gcc
CPPFLAGS=-fPIC -Wall
LDFLAGS=-shared -Wl,-soname,test.so -lc

PARTS=libfreckle.o DotStore.o DotStorageChunk.o

all: libfreckle.so

tests: runtests

libfreckle.o: libfreckle.cpp
	$(CPP) $(CPPFLAGS) -c libfreckle.cpp

libfreckle.so: $(PARTS)
	$(CPP) $(CPPFLAGS) $(LDFLAGS) $(PARTS) -o libfreckle.so


DotStorageChunk.o: DotStorageChunk.cpp DotStorageChunk.h
	$(CPP) $(CPPFLAGS) -c DotStorageChunk.cpp

DotStore.o: DotStore.cpp DotStore.h
	$(CPP) $(CPPFLAGS) -c DotStore.cpp


#
# CPP Test Suites
#
testDotStorageChunk.cpp: testDotStorageChunk.h DotStorageChunk.cpp DotStorageChunk.h
	./cxxtestgen.pl --error-printer -o testDotStorageChunk.cpp testDotStorageChunk.h

testDotStorageChunk: testDotStorageChunk.cpp DotStorageChunk.o
	$(CPP) -I./ -o testDotStorageChunk testDotStorageChunk.cpp DotStorageChunk.o

testDotStore.cpp: testDotStore.h DotStore.cpp DotStore.h
	./cxxtestgen.pl --error-printer -o testDotStore.cpp testDotStore.h

testDotStore: testDotStore.cpp DotStore.o DotStorageChunk.o
	$(CPP) -I./ -o testDotStore testDotStore.cpp DotStore.o DotStorageChunk.o

runtests: testDotStorageChunk testDotStore
	./testDotStorageChunk
	./testDotStore



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