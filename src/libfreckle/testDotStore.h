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
// 	void testDelMany(void)
// 	{
// 		DotStore *ds=new DotStore();
// 
// 		for(int i=0; i<TEST_DOTSTORE_NUMPOINTS; i++)
// 		{
// 			ds->AddDot(i,i,i);
// 			TS_ASSERT(ds->GetNum()==i+1);
// 		}
// 
// 		// test the top point is correct then delete it
// 		for(int i=0; i<TEST_DOTSTORE_NUMPOINTS; i++)
// 		{
// 			TS_ASSERT(ds->GetDot(0)->x==i);
// 			ds->DelDot(0);
// 		}
// 
// 		// free the dotstore
// 		delete ds;
// 	}

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

	// test CountAreaMatches
	void testCountAreaMatches(void)
	{
		DotStore *ds=new DotStore();

		// one big line down the center
		for(int i=0; i<300; i++)
			ds->AddDot(i,i,300-i);

		// index
		ds->CreateIndex();
// 		ds->DumpIndex();

		// lets count up sections
		// check for diagonal mid section
		TS_ASSERT(ds->CountAreaMatches(100,100,200,200,10)==100);

		// test a top left corner with window
		TS_ASSERT(ds->CountAreaMatches(0,0,10,10,10)==10);

		// test an offset block which just touches the diagonal is truly empty
		TS_ASSERT(ds->CountAreaMatches(10,0,20,10,10)==0);

		// test an offset block which intersects the bottom left is correct
		TS_ASSERT(ds->CountAreaMatches(10,1,20,11,10)==1);
		TS_ASSERT(ds->CountAreaMatches(10,5,20,15,10)==5);

		// test an offset block which just intersects thr top right is correct
		TS_ASSERT(ds->CountAreaMatches(10,20,25,35,10)==5);

		// test the bottom right corner
		TS_ASSERT(ds->CountAreaMatches(250,250,300,300,10)==50);

		ds->DestroyIndex();

		delete ds;
	}

	// test CalculateAverageGrid
	void testCalculateAverageGrid(void)
	{
		DotStore *ds=new DotStore();

		// one big line down the center
		for(int i=0; i<1000; i++)
			ds->AddDot(i,i,1000-i);

		// index
		ds->CreateIndex();
		
		int *storage=ds->CalculateAverageGrid(1000,1000,100,10);
		TS_ASSERT(storage);

		for(int y=0; y<100; y++)
			for(int x=0; x<100; x++)
				if(x==y)
					TS_ASSERT(storage[y*100+x]==10);
		
		TS_ASSERT(ds->GetMinimumGridValue(100,100,storage)==0);
		TS_ASSERT(ds->GetMaximumGridValue(100,100,storage)==10);

		ds->DestroyIndex();

		delete ds;
	}

	// testGridToString
	void testGridToString(void)
	{
		DotStore *ds=new DotStore();

		// one big line down the center
		for(int i=0; i<1000; i++)
			ds->AddDot(i,i,1000-i);

		// index
		ds->CreateIndex();
		
		int *storage=ds->CalculateAverageGrid(1000,1000,100,10);
		TS_ASSERT(storage);

		unsigned char *string=ds->GridToString();

		for(int y=0; y<100; y++)
			for(int x=0; x<100; x++)
				TS_ASSERT( string[y*100+x]==(x==y?255:0) );
		
		ds->DestroyIndex();

		delete ds;
	}

	// this data was causing a segfault. So we recreate it here. to make sure theses no segfault
	void testRealProblemData(void)
	{
		DotStore *ds=new DotStore();

		int data[226][3]={{79,4,5},{93,7,5},{87,7,5},{64,8,7},{59,8,6},{57,8,8},{65,9,6},{60,9,5},{58,9,7},{50,9,5},{64,10,6},{59,10,6},{57,10,6},{65,11,5},{60,11,5},{58,11,5},{50,11,5},{63,26,5},{61,27,10},{51,27,6},{62,28,9},{52,28,5},{63,29,8},{64,30,7},{59,30,6},{57,30,8},{65,31,6},{60,31,5},{58,31,7},{50,31,5},{64,32,6},{59,32,6},{57,32,8},{65,33,5},{60,33,5},{58,33,7},{50,33,5},{64,34,5},{59,34,6},{57,34,5},{41,40,6},{36,42,5},{91,43,5},{85,43,5},{41,45,10},{42,46,9},{63,49,5},{94,50,5},{88,50,5},{66,50,5},{47,51,5},{4,51,7},{74,52,6},{70,52,6},{25,52,6},{21,52,6},{17,52,6},{13,52,6},{9,52,6},{5,52,6},{75,53,5},{71,53,5},{22,53,5},{18,53,5},{14,53,5},{10,53,5},{6,53,5},{55,54,8},{56,55,7},{64,56,6},{59,56,11},{57,56,6},{65,57,5},{60,57,10},{58,57,5},{50,57,7},{61,58,9},{51,58,6},{62,59,8},{52,59,5},{63,60,7},{64,61,6},{59,61,7},{57,61,7},{65,62,5},{60,62,6},{58,62,6},{50,62,6},{64,63,5},{59,63,5},{57,63,6},{1,64,5},{54,70,5},{76,71,7},{72,71,5},{23,71,5},{19,71,5},{15,71,5},{11,71,5},{7,71,5},{77,72,6},{96,73,7},{97,74,6},{98,75,5},{54,78,5},{76,79,5},{72,79,9},{23,79,7},{19,79,11},{15,79,10},{11,79,10},{7,79,10},{80,80,5},{73,80,8},{24,80,6},{20,80,10},{16,80,9},{12,80,9},{8,80,9},{74,81,7},{70,81,8},{25,81,5},{21,81,9},{17,81,8},{13,81,8},{9,81,8},{5,81,8},{75,82,6},{71,82,7},{22,82,8},{18,82,7},{14,82,7},{10,82,7},{6,82,7},{76,83,5},{72,83,6},{23,83,7},{19,83,6},{15,83,6},{11,83,6},{7,83,6},{33,85,8},{34,86,7},{35,87,6},{36,88,5},{35,89,5},{53,99,6},{54,100,5},{76,101,11},{72,101,6},{23,101,8},{19,101,6},{15,101,6},{11,101,6},{7,101,6},{80,102,7},{73,102,5},{24,102,7},{20,102,5},{16,102,5},{12,102,5},{8,102,5},{81,103,6},{48,103,5},{45,103,8},{82,104,5},{46,104,7},{47,105,6},{4,105,5},{81,106,7},{48,106,5},{45,106,7},{82,107,6},{46,107,6},{95,108,6},{89,108,5},{83,108,5},{78,109,6},{79,110,5},{77,111,6},{96,112,10},{97,113,9},{98,114,8},{99,115,7},{40,115,5},{100,116,6},{101,117,5},{80,118,8},{73,118,5},{24,118,6},{20,118,5},{16,118,5},{12,118,5},{8,118,5},{81,119,7},{48,119,5},{45,119,7},{82,120,6},{46,120,6},{95,121,6},{89,121,5},{83,121,5},{78,122,7},{79,123,6},{80,124,5},{73,124,7},{24,124,6},{20,124,7},{16,124,7},{12,124,7},{8,124,7},{74,125,6},{70,125,6},{25,125,5},{21,125,6},{17,125,6},{13,125,6},{9,125,6},{5,125,6},{75,126,5},{71,126,5},{22,126,5},{18,126,5},{14,126,5},{10,126,5},{6,126,5},};

		for(int i=0; i<226; i++)
		{
			ds->AddDot(data[i][0],data[i][1],data[i][2]);
		}

		// now create our index
		ds->CreateIndex();

		int *storage=ds->CalculateAverageGrid(106,131,10,4);

		ds->DestroyIndex();
	
		delete ds;

	}
};


