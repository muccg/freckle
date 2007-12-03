#include "DotStorageChunk.h"
#include <malloc.h>			// gives us NULL
#include <assert.h>

DotStorageChunk::DotStorageChunk()
{
	num=0;
	next=NULL;
	prev=NULL;
}

DotStorageChunk::~DotStorageChunk()
{

}

Dot *DotStorageChunk::GetDot(int index)
{
	assert(index>=0);
	assert(index<DOTSTORAGECHUNKSIZE);
	return &dots[index];
}

void DotStorageChunk::AddDot(int x, int y, int length)
{
	assert(num<DOTSTORAGECHUNKSIZE);		// make sure we are not full
	dots[num].x=x;
	dots[num].y=y;
	dots[num].length=length;
	num++;
}

void DotStorageChunk::DelDot(int index)
{
	assert(index<num);
	// move them all down one
	for(int i=index; i<num; i++)
	{
		dots[i].x=dots[i+1].x;
		dots[i].y=dots[i+1].y;
		dots[i].length=dots[i+1].length;
	}
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


