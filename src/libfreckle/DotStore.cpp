#include "DotStore.h"
#include <assert.h>
#include <memory.h>
#include <stdio.h>

//construct
DotStore::DotStore()
{
	head=NULL;
	tail=NULL;

	numchunks=0;
	numdots=0;
}

//destruct
DotStore::~DotStore()
{
	// free all allocated chunks
	DotStorageChunk *next=NULL;
	for(DotStorageChunk *chunk=head; chunk; chunk=next)
	{
		next=chunk->GetNext();			// gotta save this 
		delete chunk;
	}
}

void DotStore::AddDotStorageChunk()
{
	if(!head && !tail)
	{
		//first chunk
		head=tail=new DotStorageChunk();
		assert(head);
		assert(tail);
		numchunks++;
	}
	else
	{
		//subsequenct chunk
		DotStorageChunk *newchunk=new DotStorageChunk();
		assert(newchunk);
		
		tail->LinkAfter(newchunk);
		tail=newchunk;
		numchunks++;
	}
}

// This will collapse the dot storage chunks into the minimum size
// so it will align all the dots into the lowest chunks and free any unused chunks on top
void DotStore::CollapseDotStorageChunks()
{
	assert(0);		// Not Implemented
}


// AddDot(x,y,length)
void DotStore::AddDot(int x, int y, int length)
{
	// if we have no chunks
	if(!head && !tail)
		AddDotStorageChunk();

	//find the first non empty chunk
	DotStorageChunk *chunk=FindFirstNonEmptyChunk();
	
	if(chunk)
	{
		// chunk with space found
		chunk->AddDot(x,y,length);
		numdots++;
	}
	else
	{
		// no available chunks. Add chunk and repeat
		AddDotStorageChunk();
		return AddDot(x,y,length);		//recurse
	}
}

Dot *DotStore::GetDot(int index)
{
	assert(index>=0);
	if(!head || !tail || !numdots || index>numdots-1)
		return NULL;			//no dots here, or index too high or too low

	DotStorageChunk *chunk=head;
	
	// return the pointer to dot with index 'index'
	while(index>=chunk->GetNum())
	{
		index-=chunk->GetNum();
		chunk=chunk->GetNext();
		assert(chunk);			// if we run out of chunks something horrid has happened. Because above we tested if index was too high
	}

	//so chunk now contains the relevant dot. return it
	return chunk->GetDot(index);
}

void DotStore::DelDot(int index)
{
	assert(index>=0);
	assert(head && tail && numdots && index<numdots);		//no dots here, or index too high or too low

	DotStorageChunk *chunk=head;
	
	// find the dot with index 'index'
	while(index>=chunk->GetNum())
	{
		index-=chunk->GetNum();
		chunk=chunk->GetNext();
		assert(chunk);			// if we run out of chunks something horrid has happened. Because above we tested if index was too high
	}

	//so chunk now contains the relevant dot. lets delete it
	chunk->DelDot(index);
	numdots--;
}

void DotStore::Dump()
{
	printf("i\tx\ty\tlen\n=\t=\t=\t===\n");
	for(int i=0; i<numdots; i++)
		printf("%d\t%d\t%d\t%d\n",i,GetDot(i)->x,GetDot(i)->y,GetDot(i)->length);
}



