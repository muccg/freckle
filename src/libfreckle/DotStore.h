#ifndef _DOTSTORE_H_
#define _DOTSTORE_H_

#include "Dot.h"
#include "DotStorageChunk.h"
#include "QuadTree.h"

// Our dot storage class
class DotStore
{
private:
	DotStorageChunk *head, *tail;
	int numchunks;
	int numdots;

	int maxx, maxy;
	
	void AddDotStorageChunk();
	void CollapseDotStorageChunks();

public:
	DotStore();
	~DotStore();
	
	void AddDot(int x, int y, int len);
	Dot *GetDot(int index);
	void DelDot(int index);
	
	//! \brief empty the entire store of all its dots
	void Empty();

	/**
	** \brief Turns the dotstore into a buffer for saving to disk/ram/whatever
	*/
	int *ToBuffer();
	
	//!
	//! \brief Fills the dotstore by decoding the passed in buffer
	void FromBuffer(int *buffer);

	//! \brief returns the size in bytes of the buffer that is passed in
	int BufferSize(int *buffer);

	//! \brief filter out any dots that are less than a particular length
	DotStore *Filter(int minlength);
	
	//! \brief interpolate long matches into many small matches
	void Interpolate(int window);

	inline int GetMaxX() const
	{
		return maxx;
	}
	
	inline int GetMaxY() const
	{
		return maxy;
	}

	inline int SetMaxX(int m)
	{
		maxx=m;
	}
	
	inline int SetMaxY(int m)
	{
		return maxy=m;
	}

	
	
	inline int GetDotX(int index)
	{
		return GetDot(index)->x;
	}

	inline int GetDotY(int index)
	{
		return GetDot(index)->y;
	}

	inline int GetDotLength(int index)
	{
		return GetDot(index)->length;
	}

	inline int GetNum()
	{
		return numdots;
	}

	// find first non empty chunk
	inline DotStorageChunk *FindFirstNonEmptyChunk()
	{
		assert(head);
		for(DotStorageChunk *chunk=head; chunk; chunk=chunk->GetNext())
			if(!chunk->IsFull())
				return chunk;

		// none found
		return NULL;
	}

	// Dump out the contents. For debug mainly.
	void Dump();

/*
**
** Sorted Indexing Functions
** =========================
** the following functions are to create a sorted rapid access index of all the dots once they have been calculated.
** 
** It works like this
** for each y value stored sorted in an array there is a pointer
** following this pointer gives us an array of x values sorted in order
** with each value is the length of the match
*/
private:
	QuadTree *index;
	int pixwidth;
	int pixheight;
	int *averagearray;

public:
	void CreateIndex();
	void DestroyIndex();
	void DumpIndex();

	// index access functions
	Dot *GetIndexDot(int x, int y);

	//
	// \brief sum the amount of dots within the passed in window 
	// 
	// uses an efficient algorithm and the index must have been created
	//
	int CountAreaMatches(double x1, double x2, double y1, double y2, int window);

};

#endif
