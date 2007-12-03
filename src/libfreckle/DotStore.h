#ifndef _DOTSTORE_H_
#define _DOTSTORE_H_

#include "DotStorageChunk.h"

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
};

#endif
