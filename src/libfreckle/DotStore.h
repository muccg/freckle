#ifndef _DOTSTORE_H_
#define _DOTSTORE_H_

#include "DotStorageChunk.h"
#include "LinkedListVal.h"

/* for indexing */
struct structIndexXNode
{
	int x;							//our x position 
	Dot *dot;						//our dot
};
typedef struct structIndexXNode IndexXNode;

struct structIndexYNode
{
	int y;							// our y position
	LinkedListVal<struct structIndexXNode *> *child;		// for x nodes
};
typedef struct structIndexYNode IndexYNode;




// Our dot storage class
class DotStore
{
private:
	DotStorageChunk *head, *tail;
	int numchunks;
	int numdots;
	
	void AddDotStorageChunk();
	void CollapseDotStorageChunks();

public:
	DotStore();
	~DotStore();
	
	void AddDot(int x, int y, int len);
	Dot *GetDot(int index);
	void DelDot(int index);

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
	LinkedListVal<IndexYNode *> *index;

public:
	void CreateIndex();
	void IndexDot( Dot *dot );
	void DestroyIndex();
	void DumpIndex();

	// index access functions
	Dot *GetIndexDot(int x, int y);
};

#endif
