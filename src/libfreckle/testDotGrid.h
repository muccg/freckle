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
		DotGrid *ds=new DotGrid();
		TS_ASSERT(ds);
		delete ds;
	}

};


