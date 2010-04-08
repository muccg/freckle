CPP=g++
CC=gcc

# debug
# CPPFLAGS=-fPIC -Wall -g -pg
# LDFLAGS=-shared -Wl,-soname,test.so -g -pg -lc 

#CPPFLAGS=-fPIC -Wall -march=opteron
ARCH = opteron

# not debug
CPPFLAGS=-fPIC -march=$(ARCH) -Wall -m64
LDFLAGS=-shared -Wl,-soname,test.so -march=$(ARCH) -Wall -m64

PARTS=libfreckle.o DotStore.o DotStorageChunk.o DotGrid.o QuadTreeNode.o QuadTree.o

INSTALLVERSION=0.2
TARGET=/usr/local/lib
HEADERTARGET=/usr/local/include

all: libfreckle.so

install: libfreckle.so
	-cp libfreckle.so $(TARGET)/libfreckle.so.$(INSTALLVERSION)
	-rm $(TARGET)/libfreckle.so
	-ln -s $(TARGET)/libfreckle.so.$(INSTALLVERSION) $(TARGET)/libfreckle.so
	-cp libfreckle.h $(HEADERTARGET)/
	-ldconfig

tests: runtests

libfreckle.o: libfreckle.cpp
	$(CPP) $(CPPFLAGS) -c libfreckle.cpp

libfreckle.so: $(PARTS)
	$(CPP) $(CPPFLAGS) $(LDFLAGS) $(PARTS) -o libfreckle.so


DotStorageChunk.o: DotStorageChunk.cpp DotStorageChunk.h
	$(CPP) $(CPPFLAGS) -c DotStorageChunk.cpp

DotStore.o: DotStore.cpp DotStore.h
	$(CPP) $(CPPFLAGS) -c DotStore.cpp

DotGrid.o: DotGrid.cpp DotGrid.h
	$(CPP) $(CPPFLAGS) -c DotGrid.cpp

LinkedList.o: LinkedList.cpp LinkedList.h
	$(CPP) $(CPPFLAGS) -c LinkedList.cpp

QuadTreeNode.o: QuadTreeNode.cpp QuadTreeNode.h
	$(CPP) $(CPPFLAGS) -c QuadTreeNode.cpp

QuadTree.o: QuadTree.cpp QuadTree.h
	$(CPP) $(CPPFLAGS) -c QuadTree.cpp






#
# CPP Test Suites
#
testDotStorageChunk.cpp: testDotStorageChunk.h DotStorageChunk.cpp DotStorageChunk.h
	./cxxtestgen.pl --error-printer -o testDotStorageChunk.cpp testDotStorageChunk.h

testDotStorageChunk: testDotStorageChunk.cpp DotStorageChunk.o
	$(CPP) $(CPPFLAGS) -I./ -o testDotStorageChunk testDotStorageChunk.cpp DotStorageChunk.o
	./testDotStorageChunk

runtestDotStorageChunk: testDotStorageChunk
	./testDotStorageChunk

testDotStore.cpp: testDotStore.h DotStore.cpp DotStore.h 
	./cxxtestgen.pl --error-printer -o testDotStore.cpp testDotStore.h

testDotStore: testDotStore.cpp DotStore.o DotStorageChunk.o QuadTree.o QuadTreeNode.o
	$(CPP) $(CPPFLAGS) -I./ -o testDotStore testDotStore.cpp DotStore.o DotStorageChunk.o QuadTree.o QuadTreeNode.o
# 	./testDotStore

testDotGrid.cpp: testDotGrid.h DotGrid.cpp DotGrid.h
	./cxxtestgen.pl --error-printer -o testDotGrid.cpp testDotGrid.h

testDotGrid: testDotGrid.cpp DotGrid.o DotStore.o DotStorageChunk.o QuadTree.o QuadTreeNode.o
	$(CPP) $(CPPFLAGS) -I./ -o testDotGrid testDotGrid.cpp DotGrid.o DotStore.o DotStorageChunk.o QuadTree.o QuadTreeNode.o
	./testDotGrid

testQuadTreeNode.cpp: testQuadTreeNode.h QuadTreeNode.cpp QuadTreeNode.h
	./cxxtestgen.pl --error-printer -o testQuadTreeNode.cpp testQuadTreeNode.h

testQuadTreeNode: testQuadTreeNode.cpp QuadTreeNode.o QuadTree.o DotStore.o DotStorageChunk.o
	$(CPP) $(CPPFLAGS) -I./ -o testQuadTreeNode testQuadTreeNode.cpp QuadTreeNode.o DotStore.o DotStorageChunk.o QuadTree.o
	./testQuadTreeNode

runtestQuadTreeNode: testQuadTreeNode
	./testQuadTreeNode

testQuadTree.cpp: testQuadTree.h QuadTree.cpp QuadTree.h
	./cxxtestgen.pl --error-printer -o testQuadTree.cpp testQuadTree.h

testQuadTree: testQuadTree.cpp QuadTree.o DotStore.o DotStorageChunk.o QuadTreeNode.o
	$(CPP) $(CPPFLAGS) -I./ -o testQuadTree testQuadTree.cpp QuadTree.o QuadTreeNode.o DotStore.o DotStorageChunk.o

runtestQuadTree: testQuadTree
	./testQuadTree

runtests: testDotStorageChunk testDotStore testDotGrid testQuadTreeNode
	./testDotStorageChunk
	./testDotStore
	./testDotGrid
	./testQuadTreeNode



#
# Cleans
#
cleantests:
	-rm testDotStorageChunk.cpp testDotStorageChunk testDotStore.cpp testDotStore testDotGrid.cpp testDotGrid


clean: cleantests
	-rm *.so
	-rm *.o

#gcc -fPIC -g -c -Wall a.c
#gcc -fPIC -g -c -Wall b.c
#gcc -shared -Wl,-soname,libmystuff.so.1 \
#    -o libmystuff.so.1.0.1 a.o b.o -lc
