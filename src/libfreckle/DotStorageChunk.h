#ifndef _DOTSTORAGECHUNK_H_
#define _DOTSTORAGECHUNK_H_

#include <malloc.h>			//give us NULL
#include <assert.h>
#include "Dot.h"

#define DOTSTORAGECHUNKSIZE	8192

class DotStorageChunk
{
private:
	Dot	dots[DOTSTORAGECHUNKSIZE];
	
	// this is not the most efficient way RAM wise but it will do for now
	bool	inuse[DOTSTORAGECHUNKSIZE];

	int	num;

	// doubly linked list of storage chunks
	DotStorageChunk		*prev, *next;

public:
	DotStorageChunk();
	~DotStorageChunk();

	Dot	*GetDot(int num);
	void	AddDot(int x, int y, int length);
	void	DelDot(int index);

	void LinkBefore(DotStorageChunk *insert);
	void LinkAfter(DotStorageChunk *insert);

	inline int GetDotX(int num)
	{
		return GetDot(num)->x;
	}

	inline int GetDotY(int num)
	{
		return GetDot(num)->y;
	}

	inline int GetDotLength(int num)
	{
		return GetDot(num)->length;
	}

	inline void	SetNext(DotStorageChunk *chunk)
	{
		next=chunk;
	}

	inline void	SetPrev(DotStorageChunk *chunk)
	{
		prev=chunk;
	}

	inline DotStorageChunk *GetNext()
	{
		return next;
	}

	inline DotStorageChunk *GetPrev()
	{
		return prev;
	}

	inline bool	IsFull()
	{
		return num==DOTSTORAGECHUNKSIZE;
	}

	inline bool	IsEmpty()
	{
		return num==0;
	}

	inline int	GetNum()
	{
		return num;
	}
};

#endif
