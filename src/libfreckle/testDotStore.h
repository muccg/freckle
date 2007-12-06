#include <cxxtest/TestSuite.h>

#include "DotStore.h"

class MyTestSuite : public CxxTest::TestSuite
{
public:
	//
	// Test construction and freeing of the class object
	//
	void testConstructFree(void)
	{
		DotStore *ds=new DotStore();
		TS_ASSERT(ds);
		delete ds;
	}

	// test adding tens of thousands of points
	void testAddMany(void)
	{
		#ifndef TEST_DOTSTORE_NUMPOINTS
		#define TEST_DOTSTORE_NUMPOINTS	100000
		#endif

		DotStore *ds=new DotStore();

		for(int i=0; i<TEST_DOTSTORE_NUMPOINTS; i++)
		{
			ds->AddDot(i,i,i);
			TS_ASSERT(ds->GetNum()==i+1);
		}

		// test they are all there and correct
		for(int i=0; i<TEST_DOTSTORE_NUMPOINTS; i++)
			TS_ASSERT(ds->GetDot(i)->x==i);

		// free the dotstore
		delete ds;
	}

	// test adding tens of thousands and then deleting them
	void testDelMany(void)
	{
		DotStore *ds=new DotStore();

		for(int i=0; i<TEST_DOTSTORE_NUMPOINTS; i++)
		{
			ds->AddDot(i,i,i);
			TS_ASSERT(ds->GetNum()==i+1);
		}

		// test the top point is correct then delete it
		for(int i=0; i<TEST_DOTSTORE_NUMPOINTS; i++)
		{
			TS_ASSERT(ds->GetDot(0)->x==i);
			ds->DelDot(0);
		}

		// free the dotstore
		delete ds;
	}

	// test indexing
	void testIndexing(void)
	{
		DotStore *ds=new DotStore();

		for(int i=0; i<TEST_DOTSTORE_NUMPOINTS; i++)
			ds->AddDot(i,int(i/100),i);

		ds->CreateIndex();
		
		//lets try and access the indexing
		TS_ASSERT(ds->GetIndexDot(20152,201)->x==20152);
		TS_ASSERT(ds->GetIndexDot(20152,201)->y==201);
		TS_ASSERT(ds->GetIndexDot(20152,201)->length==20152);

		// test some not found points
		TS_ASSERT(ds->GetIndexDot(201,201)==NULL);
		TS_ASSERT(ds->GetIndexDot(20152,20152)==NULL);
		
		// test destroying the indexing
		ds->DestroyIndex();
		
	}
};


