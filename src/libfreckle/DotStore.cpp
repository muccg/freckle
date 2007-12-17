#include "DotStore.h"
#include <assert.h>
#include <memory.h>
#include <stdio.h>
#include <math.h>

//construct
DotStore::DotStore()
{
	head=NULL;
	tail=NULL;

	numchunks=0;
	numdots=0;

	maxx=maxy=0;

	index=NULL;

	averagearray=NULL;
	pixwidth=0;
	pixheight=0;
}

//destruct
DotStore::~DotStore()
{
	// delete the index if its there
	if(index)
		DestroyIndex();

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
	// keep track of maximums
	if(x>maxx)
		maxx=x;
	if(y>maxy)
		maxy=y;

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

// indexes the dots by y and then by x
void DotStore::CreateIndex()
{
	// if there is an old index destroy it
	if(index)
		DestroyIndex();

	// create a new index
	index=new QuadTree(0,0,maxx,maxy);

	// go through every dot
	for(int i=0; i<numdots; i++)
		index->AddDot(GetDot(i));
}

void DotStore::DestroyIndex()
{
	// nuke the index
	assert(index);	

	delete index;

	index=NULL;
}


//Dump the index to stdout. Mainly for DEBUG
void DotStore::DumpIndex()
{
	assert(index);
	index->Dump();
}


// use the index to quickly find a dot
Dot *DotStore::GetIndexDot(int x,int y)
{
	assert(index);				//we must be indexed
	
	LinkedListVal<Dot *> *result=index->SpatialQuery(x,y,x,y);
	assert(result->Length()==1 || result->Length()==0);

	if(result->Length()==0)
		return NULL;

	Dot *dot=result->Pop();

	delete result;

	// found it
	return dot;
}

int DotStore::CountAreaMatches(double x1, double y1, double x2, double y2, int window)
{
	assert(index);				// we must be indexed
	//printf("%f,%f\n",x2-x1,y2-y1);
	//assert(x2-x1 == y2-y1);			// we must be square. TODO: add support for rectangular area

	double dwindow=(double)window;

	int count=0;

	// do a spatial query on the index and get a list of matching dots. TODO deal with fractions properly by testing the dots to make sure they should *really* be included
	LinkedListVal<Dot *> *result=index->SpatialQuery((int)floor(x1-dwindow),(int)floor(y1-dwindow),(int)ceil(x2),(int)ceil(y2));

	for(LinkedListVal<Dot *>::Iterator i(*result); !i.Done(); i++)
	{
		if ( (*i)->x >= x1-dwindow && (*i)->y >= y1-dwindow && (*i)->x < x2 && (*i)->y < y2 )
		{ 
			// applicable x and y
			double x=(double)(*i)->x+0.5;
			double y=(double)(*i)->y+0.5;
			double length=(double)(*i)->length;
			double protrude=length;			//how much protrudes into this calculation square
			// which zone are we in?
			// 1. inside the window
			if( x>=x1 && x<x2 && y>=y1 && y<y2)
			{
				//we are inside the window
				// truncate protrude if we extend outside to the bottom or the right of the window
				if(x+protrude > x2)
					protrude=x2-x;
				if(y+protrude > y2)
					protrude=y2-y;
				
	
				//scan down and to the right to see when our next match point comes up or until our length is exhausted. add one for each point
				if(protrude)
				{
					int xp=(int)x;
					int yp=(int)y;
					do
					{ 
						count++;
						protrude-=1.0;
					} while( (!GetIndexDot(++xp,++yp)) && protrude>=1.0);
				}
			}
			// 2. The parallelogram above the areaand including right on the line
			else if( x>=(y-y1+x1) && x<(y-y1+x2) )
			{
				//we are above the window
				if(length>y1-y)
				{
					//we extend into the window
					protrude=length-(y1-y);
					
					
					// check if we extend out of the window to the right too
					double sigma=length-(x2-x);
					if(sigma>0)
						protrude-=sigma;
			
					// TODO: Check if we extend out of the bottom (non square window)
	
					int xp=(int)x;
					int yp=(int)y;
					do
					{
						if(yp>=y1)
						{
							count++;
							protrude-=1.0;
						}
					} while((!GetIndexDot(++xp,++yp)) && protrude >=1.0);
				}
			}
			// 3. The parallelogram to the left of the area
			else if( y>(x-x1+y1) && y<(x-x1+y2) )
			{
				//we are to the left of the window
				if(length>x1-x)
				{
					protrude=length-(x1-x);
					
					// check if we extend out of the window to the bottom too
					double sigma=length-(y2-y);
					if(sigma>0)
						protrude-=sigma;
						
					// TODO: Check if we extend out of the window to the right, too (non square window)
					
					int xp=(int)x;
					int yp=(int)y;
					do
					{
						if(xp>=x1)
						{
							count++;
							protrude-=1.0;
						}
					} while((!GetIndexDot(++xp,++yp)) && protrude >= 1.0);
				}
			}
		}
	}

	delete result;

	return count;
}
