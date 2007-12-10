#include <cxxtest/TestSuite.h>

#include "DotGrid.h"

class MyTestSuite : public CxxTest::TestSuite
{
public:
	//
	// Test construction and freeing of the class object
	//
	void testConstructFree(void)
	{
		DotGrid *dg=new DotGrid();
		TS_ASSERT(dg);
		delete dg;
	}

	void testCreateDestroy(void)
	{
		DotGrid *dg=new DotGrid();
		TS_ASSERT(dg);

		dg->Create(101,102);
		TS_ASSERT(dg->GetWidth()==101);
		TS_ASSERT(dg->GetHeight()==102);
		TS_ASSERT(dg->GetSize()==102*101);
		dg->Destroy();

		delete dg;
	}

};


