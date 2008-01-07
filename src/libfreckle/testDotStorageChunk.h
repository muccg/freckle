#include <cxxtest/TestSuite.h>

#include "DotStorageChunk.h"

class MyTestSuite : public CxxTest::TestSuite
{
public:
	//
	// Test construction and freeing of the class object
	//
	void testConstructFree(void)
	{
		DotStorageChunk *dsc=new DotStorageChunk();
		TS_ASSERT(dsc);
		delete dsc;
	}

	// Test nexts and previouses
	void testLink(void)
	{
		DotStorageChunk *chunk1=new DotStorageChunk();
		DotStorageChunk *chunk2=new DotStorageChunk();
		TS_ASSERT(chunk1);
		TS_ASSERT(chunk2);

		// They should begin being NULL
		TS_ASSERT(chunk1->GetPrev()==NULL);
		TS_ASSERT(chunk2->GetPrev()==NULL);
		TS_ASSERT(chunk1->GetNext()==NULL);
		TS_ASSERT(chunk2->GetNext()==NULL);

		// now test the setting and getting of the chunks
		chunk1->SetNext(chunk2);
		TS_ASSERT(chunk1->GetNext()==chunk2);
		chunk1->SetPrev(chunk2);
		TS_ASSERT(chunk1->GetPrev()==chunk2);
		chunk1->SetNext(NULL);
		TS_ASSERT(chunk1->GetNext()==NULL);
		chunk1->SetPrev(NULL);
		TS_ASSERT(chunk1->GetPrev()==NULL);

		// Lets try the linking code
		chunk1->LinkAfter(chunk2);
		TS_ASSERT(chunk1->GetNext()==chunk2);
		TS_ASSERT(chunk2->GetPrev()==chunk1);

		// delink
		chunk1->LinkAfter(NULL);
		TS_ASSERT(chunk1->GetNext()==NULL);
		TS_ASSERT(chunk2->GetPrev()==NULL);

		// link before
		chunk1->LinkBefore(chunk2);
		TS_ASSERT(chunk1->GetPrev()==chunk2);
		TS_ASSERT(chunk2->GetNext()==chunk1);

		// delink
		chunk1->LinkBefore(NULL);
		TS_ASSERT(chunk1->GetPrev()==NULL);
		TS_ASSERT(chunk2->GetNext()==NULL);
		
		// free resources
		delete chunk1;
		delete chunk2;
	}
	
	void testAddDot()
	{
		DotStorageChunk *chunk=new DotStorageChunk();
		TS_ASSERT(chunk);

		int count=0;
		TS_ASSERT(chunk->IsEmpty());

		while(count<DOTSTORAGECHUNKSIZE)
		{
			count++;
			chunk->AddDot(10,-40,100);
			TS_ASSERT(chunk->GetNum()==count);
		}
		
		TS_ASSERT(chunk->IsFull());
		
		delete chunk;
	}

	void testDelDot()
	{
		DotStorageChunk *chunk=new DotStorageChunk();
		TS_ASSERT(chunk);

		//first lets fill it up
		int xyz=0;
		while(!chunk->IsFull())
		{
			chunk->AddDot(xyz,xyz,5);
			xyz++;
		}
	
		//we are now full
		int num=chunk->GetNum();
		TS_ASSERT(num==DOTSTORAGECHUNKSIZE);

		//delete from the head one by one
		Dot *dp=chunk->GetDot(0);		// our head dot. We will watch this as we destroy the list
		for(int i=0; i<num; i++)
		{
			TS_ASSERT(dp->x==i);
			TS_ASSERT(dp->y==i);
			TS_ASSERT(dp->length==5);

			//delete the head
			chunk->DelDot(0);
			dp=chunk->GetDot(0);
		}
		TS_ASSERT(chunk->IsEmpty());

		delete chunk;
	}

	void testAddition( void )
	{
		TS_ASSERT( 1 + 1 > 1 );
		TS_ASSERT_EQUALS( 1 + 1, 2 );
	}
};

