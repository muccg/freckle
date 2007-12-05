#ifndef _DOTSTORE_H_
#define _DOTSTORE_H_

#include "DotStorageChunk.h"

/* for indexing */
struct structIndexNode
{
	struct structIndexNode *next;					// next node in this node array
	int index;							// our value
	union
	{
		struct structIndexNode *child;				// for nodes
		int length;						// for leaves
	} value;
};
typedef struct structIndexNode IndexNode;

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
	IndexNode *yhead,*ytail;

public:
	void CreateIndex();
	
};

#endif
