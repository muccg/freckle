#include <cxxtest/TestSuite.h>

#include "QuadTreeNode.h"

class MyTestSuite : public CxxTest::TestSuite
{
public:
	//
	// Test construction and freeing of the class object
	//
	void testConstructFree(void)
	{
		QuadTreeNode *qtn=new QuadTreeNode();
		TS_ASSERT(qtn);
		delete qtn;
	}

	#define NUMX 100
	#define NUMY 100
	void testAddDots(void)
	{
		Dot dot[NUMX*NUMY];
		for(int x=0; x<NUMX; x++)
			for(int y=0; y<NUMY; y++)
			{
				dot[y*NUMX+x].x=x;
				dot[y*NUMX+x].y=y;
				dot[y*NUMX+x].length=1;
			}

		QuadTreeNode *qtn=new QuadTreeNode(0,0,NUMX,NUMY);
		TS_ASSERT(qtn);
		
		for(int i=0; i<NUMX*NUMY; i++)
			qtn->AddDot( &dot[i] );

		// do some simple spatial queries to make sure we have what we expect
		LinkedListVal<Dot *> *list=new LinkedListVal<Dot *>;
		qtn->SpatialQueryRecurse(list,10,10,19,19);
		TS_ASSERT(list->Length()==100);
		delete list;

		list=new LinkedListVal<Dot *>;
		qtn->SpatialQueryRecurse(list,15,34,89,59);
		TS_ASSERT(list->Length()==(89-15+1)*(59-34+1));
		delete list;

		list=new LinkedListVal<Dot *>;
		qtn->SpatialQueryRecurse(list,0,0,100,100);
		TS_ASSERT(list->Length()==10000);
		delete list;

		delete qtn;
	}

	void testDelDots1(void)
	{
		Dot dot[NUMX*NUMY];
		for(int x=0; x<NUMX; x++)
			for(int y=0; y<NUMY; y++)
			{
				dot[y*NUMX+x].x=x;
				dot[y*NUMX+x].y=y;
				dot[y*NUMX+x].length=1;
			}

		QuadTreeNode *qtn=new QuadTreeNode(0,0,NUMX,NUMY);
		TS_ASSERT(qtn);
		
		for(int i=0; i<NUMX*NUMY; i++)
			qtn->AddDot( &dot[i] );

		// now we try deleting a quarter of the dots
		for(int x=0; x<NUMX; x++)
			for(int y=0; y<NUMY; y++)
				if((x%2==0) && (y%2==0))
					qtn->DeleteDot( &dot[y*NUMX+x] );


		// spatial query of whole area should return half the number now
		LinkedListVal<Dot *> *list=new LinkedListVal<Dot *>;
		qtn->SpatialQueryRecurse(list,0,0,100,100);
		TS_ASSERT(list->Length()==7500);
		delete list;

		delete qtn;
	}

	void testDelDots2(void)
	{
		Dot dot[NUMX*NUMY];
		for(int x=0; x<NUMX; x++)
			for(int y=0; y<NUMY; y++)
			{
				dot[y*NUMX+x].x=x;
				dot[y*NUMX+x].y=y;
				dot[y*NUMX+x].length=1;
			}

		QuadTreeNode *qtn=new QuadTreeNode(0,0,NUMX,NUMY);
		TS_ASSERT(qtn);
		
		for(int i=0; i<NUMX*NUMY; i++)
			qtn->AddDot( &dot[i] );

		// now we try deleting every second dot
		for(int i=0; i<NUMX*NUMY;i+=2)
			qtn->DeleteDot( &dot[i] );


		// spatial query of whole area should return half the number now
		LinkedListVal<Dot *> *list=new LinkedListVal<Dot *>;
		qtn->SpatialQueryRecurse(list,0,0,100,100);
		TS_ASSERT(list->Length()==5000);
		delete list;

		delete qtn;
	}
};



