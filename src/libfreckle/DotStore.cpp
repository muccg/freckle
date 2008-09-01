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
	Empty();
}

//empty
void DotStore::Empty()
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

Dot *DotStore::GetDot(int ind)
{
	assert(ind>=0);
	if(!head || !tail || !numdots || ind>numdots-1)
		return NULL;			//no dots here, or index too high or too low

	DotStorageChunk *chunk=head;
	
	// return the pointer to dot with index 'index'
	while(ind>=chunk->GetNum())
	{
		ind-=chunk->GetNum();
		chunk=chunk->GetNext();
		assert(chunk);			// if we run out of chunks something horrid has happened. Because above we tested if index was too high
	}

	//so chunk now contains the relevant dot. return it
	return chunk->GetDot(ind);
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
// 	Dump();
// 	
	//printf("DotStore::CreateIndex() %d %d %d\n",maxx,maxy,numdots);

	// if there is an old index destroy it
	if(index)
	{
// 		printf("destroying old index\n");
		DestroyIndex();
	}

	if(numdots==0)
	{
// 		printf("numdots==0\n");
		// we are an empy dotstore.
		// this is a bit of a problem because we don't know how big the sequence is. 
		// for now, we set maxx and maxy to 1 and 1
		// this shouldn't be a problem because we have not dots, so none will be added
		// really what we should do is use the input total sequence lengths here
		maxx=1;
		maxy=1;
	}

	// try to fix conserved region crash?
	if(maxx==0)
		maxx=1;

	if(maxy==0)
		maxy=1;

	// create a new index
// 	printf("creating quadtree\n");
	index=new QuadTree(0,0,maxx,maxy);

	// go through every dot
// 	printf("adding dots: %d\n",numdots);
	for(int i=0; i<numdots; i++)
	{
		index->AddDot(GetDot(i));
	}
}

void DotStore::DestroyIndex()
{
	// nuke the index
	assert(index);	

	delete index;

	index=NULL;
}

#ifdef DEBUG
//Dump the index to stdout. Mainly for DEBUG
void DotStore::DumpIndex()
{
	assert(index);
	index->Dump();
}
#endif

// use the index to quickly find a dot
Dot *DotStore::GetIndexDot(int x,int y)
{
	assert(index);				//we must be indexed
	
	LinkedListVal<Dot *> *result=index->SpatialQuery(x,y,x,y);
	assert(result->Length()==1 || result->Length()==0);

	if(result->Length()==0)
	{
		delete result;
		return NULL;
	}

	Dot *dot=result->Pop();

	delete result;

	// found it
	return dot;
}

// use the index to quickly find the longest dot match on a particular row
Dot *DotStore::GetIndexLongestMatchingRowDot(int y)
{
	assert(index);				//we must be indexed

	printf("Doing spatial query on %ld (%d,%d,%d,%d)\n",(long)index,0,y,maxx,y);
	LinkedListVal<Dot *> *result = index->SpatialQuery(0,y,maxx,y);
	
	// find out which dot matched is the longest
	Dot *longest=NULL;
	int maxlength=0;
	int dotlength=0;
	for( LinkedListVal<Dot *>::Iterator i(*result); !i.Done(); i++)
	{
		dotlength = (*i)->length;

		if(dotlength>maxlength)
		{
			maxlength=dotlength;
			longest=(*i);
		}
	}

	delete result;
	return longest;
}

// use the index to quickly find the longest dot match on a particular column
Dot *DotStore::GetIndexLongestMatchingColumnDot(int x)
{
	assert(index);				//we must be indexed

	LinkedListVal<Dot *> *result = index->SpatialQuery(x,0,x,maxy);
	
	// find out which dot matched is the longest
	Dot *longest=NULL;
	int maxlength=0;
	int dotlength=0;
	for( LinkedListVal<Dot *>::Iterator i(*result); !i.Done(); i++)
	{
		dotlength = (*i)->length;

		if(dotlength>maxlength)
		{
			maxlength=dotlength;
			longest=(*i);
		}
	}

	delete result;
	return longest;
}

int DotStore::CountAreaMatches(double x1, double y1, double x2, double y2, int window)
{
	assert(index);				// we must be indexed
	//printf("%f,%f\n",x2-x1,y2-y1);
	//assert(x2-x1 == y2-y1);			// we must be square. TODO: add support for rectangular area

	double dwindow=(double)window;

	printf("DotStore::CountAreaMatches(): %d\n",numdots);

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

int *DotStore::ToBuffer()
{
	// buffer returned is a pointer to an array of ints
	// the first two ints are the maxx and maxy records
	// the next int is the number of records (n)
	// it is followed by n*3 ints of the actual record data
	
	// allocate the buffer
	int buffsize=GetNum()*3+3;
	int *buffer=new int[buffsize];
	memset(buffer, 0, sizeof(int)*buffsize);

	int *buffp=buffer;

	// the maxx and maxy
	*buffp++=maxx;
	*buffp++=maxy;

	// our record number
	*buffp++=GetNum();

	// now all the records
	Dot *dot;
	for(int i=0; i<GetNum(); i++)
	{
		dot=GetDot(i);
		*buffp++=dot->x;
		*buffp++=dot->y;
		*buffp++=dot->length;
	}

	// we should perfectly fill the allocated buffer
	assert(buffp == buffer+buffsize);

	return buffer;
}

void DotStore::FromBuffer(int *buffer)
{
	// empty our dots
	Empty();

	int *buffp=buffer;
	
	// maxx and maxy
	maxx=*buffp++;
	maxy=*buffp++;

	//number of records
	int num=*buffp++;

	// read each record
	for(int i=0; i<num; i++)
	{
		AddDot(buffp[0], buffp[1], buffp[2]);
		buffp+=3;
	}

	// we should have the same number of records we expected
	assert(GetNum() == num);
}

int DotStore::BufferSize(int *buffer)
{
	return buffer[2]*3+3;
}

DotStore *DotStore::Filter(int minlength)
{
	DotStore *filteredstore=new DotStore();
	filteredstore->maxx=maxx;
	filteredstore->maxy=maxy;

	// loop through every dot and if they are less that the minlength then ignore them, else add them to the new dotstore
	Dot *dot=NULL;
	for(int i=0; i<GetNum(); i++)
	{
		dot=GetDot(i);
		if(dot->length >= minlength)
			filteredstore->AddDot(dot->x, dot->y, dot->length);
	}

	return filteredstore;
}

// Any matches that are greater than length window are processed
// added to the dot store are a bunch of sub matches to step across the window
void DotStore::Interpolate(int window)
{
	DotStore *extradots=new DotStore();

	Dot *dot=NULL;

	printf("PREInterpolate\n");
	Dump();

	for(int i=0; i<GetNum(); i++)
	{
		dot=GetDot(i);

		if(dot->length > window)
		{
			// break it down, now!
			int remainder=dot->length-window;
			int xpos=dot->x+window;
			int ypos=dot->y+window;

			while(remainder>0)
			{
				extradots->AddDot(xpos,ypos,remainder);
				xpos+=window;
				ypos+=window;
				remainder-=window;
			}
		}
	}

	// now add all these extra dots into our store
	for(int i=0; i<extradots->GetNum(); i++)
	{
		dot=extradots->GetDot(i);
		AddDot(dot->x,dot->y,dot->length);
	}

	delete extradots;

	printf("POSTInterpolate\n");
	Dump();
}
