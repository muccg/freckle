#include "DotStorageChunk.h"
#include <malloc.h>			// gives us NULL
#include <assert.h>
#include <stdio.h>

DotStorageChunk::DotStorageChunk()
{
	num=0;
	
	next=NULL;
	prev=NULL;

	for(int i=0; i<DOTSTORAGECHUNKSIZE; i++)
		inuse[i]=false;
}

DotStorageChunk::~DotStorageChunk()
{

}

Dot *DotStorageChunk::GetDot(int index)
{
	assert(index>=0);
	assert(index<DOTSTORAGECHUNKSIZE);

	// this array may be spacious, so we have to count in past unused spots to find the indexed dot
	int n=0;
	int ind=0;
	while(ind<=index and n-1<DOTSTORAGECHUNKSIZE)
	{
		if(inuse[n])
		{
			ind++;
		}
		n++;
	}

	if(n-1 == DOTSTORAGECHUNKSIZE)
		return NULL;

	return &dots[n-1];
}

void DotStorageChunk::AddDot(int x, int y, int length)
{
	assert(num<DOTSTORAGECHUNKSIZE);		// make sure we are not full

	// find first available and stick it there
	int n=0;
	while(inuse[n])
		n++;
	dots[n].x=x;
	dots[n].y=y;
	dots[n].length=length;
	inuse[n]=true;					// mark as used
	num++;
}

void DotStorageChunk::DelDot(int index)
{
	assert(index<num);

	// find this index to delete
	int n=0;
	int ind=0;
	while(ind<=index)
	{
		if(inuse[n])
		{
			ind++;
		}
		n++;
	}

	assert(n-1 < DOTSTORAGECHUNKSIZE);
	inuse[n-1]=false;				// mark as unused
	
	num--;	
}

// link 'insert' in before this in any linked chain
void DotStorageChunk::LinkBefore(DotStorageChunk *insert)
{
	//if you pass NULL in it breaks the before link
	if(!insert)
	{
		if(GetPrev())
			GetPrev()->SetNext(NULL);
		SetPrev(NULL);
	}
	else
	{
		//insert has to be standing alone
		assert(insert->GetPrev()==NULL);
		assert(insert->GetNext()==NULL);

		//link it in
		insert->SetPrev(GetPrev());
		if(GetPrev())			// there was someone before (!NULL)
			GetPrev()->SetNext(insert);
		insert->SetNext(this);
		SetPrev(insert);
	}
}

// link insert in after this in any linked chain
void DotStorageChunk::LinkAfter(DotStorageChunk *insert)
{
	//if you pass in NULL it breaks the after link
	if(!insert)
	{
		if(GetNext())
			GetNext()->SetPrev(NULL);
		SetNext(NULL);
	}
	else
	{
		//insert has to be standing alone
		assert(insert->GetPrev()==NULL);
		assert(insert->GetNext()==NULL);
	
		//link it in
		insert->SetNext(GetNext());
		if(GetNext())			// there was something after (!NULL)
			GetNext()->SetPrev(insert);
		insert->SetPrev(this);
		SetNext(insert);
	}
}


