#include <cxxtest/TestSuite.h>

#include "DotStore.h"

#include <stdlib.h>

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

		for(int i=0; i<TEST_DOTSTORE_NUMPOINTS/4; i++)
		{
			ds->AddDot(i,i,i);
			TS_ASSERT(ds->GetNum()==i+1);
		}

		// test the top point is correct then delete it
		for(int i=0; i<TEST_DOTSTORE_NUMPOINTS/4; i++)
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

	// if we compare any sequence, then ask for an average of any square that has the averaging area
	// line up with a single nucleotide base pair exactly (x,y - x+1,y+1 where x and y are integers)
	// then the value returned MUST be 0 or 1
	void testAreaAverageBoundaryCondition(void)
	{
		DotStore *ds=new DotStore();

		int data[198][3]={{31,0,6},{0,0,114},{1,1,113},{2,2,112},{3,3,111},{4,4,110},{13,5,7},{5,5,109},{6,6,108},{7,7,107},{19,8,5},{8,8,106},{9,9,105},{99,10,6},{97,10,6},{47,10,5},{10,10,104},{98,11,5},{23,11,6},{11,11,103},{24,12,5},{12,12,102},{13,13,101},{5,13,7},{14,14,100},{15,15,99},{16,16,98},{17,17,97},{18,18,96},{19,19,95},{8,19,5},{20,20,94},{43,21,5},{21,21,93},{22,22,92},{98,23,6},{23,23,91},{11,23,6},{24,24,90},{12,24,5},{25,25,89},{26,26,88},{27,27,87},{30,28,5},{28,28,86},{29,29,85},{30,30,84},{28,30,5},{31,31,83},{0,31,6},{108,32,6},{32,32,82},{109,33,5},{33,33,81},{34,34,80},{35,35,79},{38,36,5},{36,36,78},{37,37,77},{38,38,76},{36,38,5},{39,39,75},{40,40,74},{41,41,73},{42,42,72},{43,43,71},{21,43,5},{44,44,70},{45,45,69},{96,46,6},{46,46,68},{99,47,6},{97,47,5},{47,47,67},{10,47,5},{100,48,5},{48,48,66},{49,49,65},{50,50,64},{51,51,63},{90,52,8},{52,52,62},{91,53,7},{53,53,61},{92,54,6},{54,54,60},{93,55,5},{55,55,59},{56,56,58},{57,57,57},{58,58,56},{59,59,55},{60,60,54},{61,61,53},{62,62,52},{65,63,5},{63,63,51},{64,64,50},{65,65,49},{63,65,5},{69,66,6},{66,66,48},{70,67,5},{67,67,47},{68,68,46},{69,69,45},{66,69,6},{70,70,44},{67,70,5},{71,71,43},{72,72,42},{73,73,41},{74,74,40},{75,75,39},{76,76,38},{77,77,37},{86,78,5},{84,78,7},{82,78,9},{80,78,11},{78,78,36},{85,79,6},{83,79,8},{81,79,10},{79,79,35},{86,80,5},{84,80,7},{82,80,9},{80,80,34},{78,80,11},{85,81,6},{83,81,8},{81,81,33},{79,81,10},{86,82,5},{84,82,7},{82,82,32},{80,82,9},{78,82,9},{85,83,6},{83,83,31},{81,83,8},{79,83,8},{86,84,5},{84,84,30},{82,84,7},{80,84,7},{78,84,7},{85,85,29},{83,85,6},{81,85,6},{79,85,6},{86,86,28},{84,86,5},{82,86,5},{80,86,5},{78,86,5},{87,87,27},{88,88,26},{89,89,25},{90,90,24},{52,90,8},{91,91,23},{53,91,7},{92,92,22},{54,92,6},{93,93,21},{55,93,5},{94,94,20},{95,95,19},{96,96,18},{46,96,6},{99,97,5},{97,97,17},{47,97,5},{10,97,6},{98,98,16},{23,98,6},{11,98,5},{99,99,15},{97,99,5},{47,99,6},{10,99,6},{100,100,14},{48,100,5},{101,101,13},{102,102,12},{107,103,5},{103,103,11},{104,104,10},{105,105,9},{106,106,8},{107,107,7},{103,107,5},{108,108,6},{32,108,6},{109,109,5},{33,109,5}};
		
		for(int i=0; i<198; i++)
		{
			ds->AddDot(data[i][0],data[i][1],data[i][2]);
		}

		// now create our index
		ds->CreateIndex();

		// now check every single pixel average is zero or one
		for(double x=0.0; x<114.0; x+=1.0)
			for(double y=0.0; y<114.0; y+=1.0)
			{
				int count=ds->CountAreaMatches(x,y,x+1.0,y+1.0,5);
				TS_ASSERT( count==0 || count==1);
			}

		ds->DestroyIndex();
	
		delete ds;

	}

	void testEmpty()
	{
		DotStore *ds=new DotStore();

		// one big line down the center
		for(int i=0; i<300; i++)
			ds->AddDot(i,i,300-i);

		// index
		ds->CreateIndex();

		// empty it
		ds->Empty();

		// check we are empty
		TS_ASSERT(ds->GetNum()==0);
		
		delete ds;
	}

	void testBufferLoadSave(void)
	{
		// create a buffer
		DotStore *ds=new DotStore();
		for(int i=0; i<258867; i++)
			ds->AddDot(i,i,258867-i);

		// save the buffer
		int *buffer=ds->ToBuffer();

		// delete the dotstore
		delete ds;

		// create a new store and load the buffer
		DotStore *ds2=new DotStore();
		ds2->FromBuffer(buffer);

		// check the points are back
		TS_ASSERT(ds2->GetNum() == 258867);
		for(int i=0; i<258867; i++)
		{
// 			printf("%d => %d,%d,%d\n",i,ds2->GetDot(i)->x,ds2->GetDot(i)->y,ds2->GetDot(i)->length);
			TS_ASSERT(ds2->GetDot(i)->x==i);
			TS_ASSERT(ds2->GetDot(i)->y==i);
			TS_ASSERT(ds2->GetDot(i)->length==258867-i);
		}

		delete ds2;
		delete buffer;
// 		delete ds;
	}

	#define RANDINT(max) (rand()%max)
	void testFilter(void)
	{
		#define NUMTESTDOTS 1000
		Dot dot[NUMTESTDOTS];
		DotStore *ds=new DotStore();
		DotStore *oldds=NULL;

		#define MAXLEN 100
		int counts[MAXLEN];
		memset(counts, 0, MAXLEN*sizeof(int));

		for(int i=0; i<NUMTESTDOTS; i++)
		{
			dot[i].x=RANDINT(1000);
			dot[i].y=RANDINT(1000);
			dot[i].length=RANDINT(MAXLEN);
			if(dot[i].length==0)
				dot[i].length=1;

			counts[dot[i].length]++;		// keep a histogram for later look up

			ds->AddDot(dot[i].x, dot[i].y, dot[i].length);
		}

		int shouldbe=NUMTESTDOTS;
		
		TS_ASSERT(ds->GetNum() == shouldbe);
		for(int j=0; j<MAXLEN; j++)
		{
			oldds=ds;
			ds=ds->Filter(j+1);			// clear out ever dot length less than OR EQUAL TO j	
			delete oldds;
			shouldbe-=counts[j];			// we should loose this many
			TS_ASSERT(ds->GetNum() == shouldbe);	// make sure thats true
		}
		
		delete ds;
	}

	// data must be in i,x,y,len columns, whitespace delimited
	#define DATASET		"segfaultdata.txt"

	void testLargeDataIndexSegfault(void)
	{
		// load our dataset
		FILE *input=NULL;
		input=fopen(DATASET, "r");
		TS_ASSERT(input!=NULL);

		// place to storethe dots
		DotStore store;

		// count the lines
		int i,x,y,len;
		int count=0;
		while(fscanf(input,"%d %d %d %d",&i,&x,&y,&len)!=EOF)
		{
			// store in dotstore
			store.AddDot(x,y,len);

			count++;
		}

		fclose(input);

// 		printf("%d dots read\n",count);

		// create index
		store.CreateIndex();
	}
};


