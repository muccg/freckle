#include <cxxtest/TestSuite.h>

#include "QuadTree.h"

#include <stdlib.h>

class MyTestSuite : public CxxTest::TestSuite
{
public:
	//
	// Test construction and freeing of the class object
	//
	void testConstructFree(void)
	{
		QuadTree *qt=new QuadTree(0,0,1000,1000);
		TS_ASSERT(qt);
		delete qt;
	}

	void testAddDots(void)
	{
		QuadTree *qt=new QuadTree(0,0,1000,1000);
		TS_ASSERT(qt);

		Dot dot[10];

		for(int i=0; i<10; i++)
		{
			dot[i].x=5*i;
			dot[i].y=4*i;
			dot[i].length=1;
			qt->AddDot(&dot[i]);
		}

		//qt->Dump();

		delete qt;
	}

	#define MANY_DOTS	100000
	void testAddManyDots(void)
	{
		QuadTree *qt=new QuadTree(0,0,350000,350000);
		TS_ASSERT(qt);

		Dot dot[MANY_DOTS];

		for(int i=0; i<MANY_DOTS; i++)
		{
			dot[i].x=rand()%350000;
			dot[i].y=rand()%350000;
			dot[i].length=rand()%100;
			qt->AddDot(&dot[i]);
		}

		//qt->Dump();

		delete qt;
	}

	void testSpatialQuery(void)
	{
		QuadTree *qt=new QuadTree(0,0,350000,350000);
		TS_ASSERT(qt);

		Dot dot[MANY_DOTS];

		for(int i=0; i<MANY_DOTS; i++)
		{
			dot[i].x=rand()%350000;
			dot[i].y=rand()%350000;
			dot[i].length=rand()%100;
			qt->AddDot(&dot[i]);
		}
	
		LinkedListVal<Dot *> *list=qt->SpatialQuery(0,0,100000,100000);
		for(LinkedListVal<Dot *>::Iterator i(*list); !i.Done(); i++)
		{	
			TS_ASSERT((*i)->x >=0 && (*i)->x <= 100000);
			TS_ASSERT((*i)->y >=0 && (*i)->y <= 100000);
		}
		delete list;
		delete qt;
	}

};

