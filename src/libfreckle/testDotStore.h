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
};


