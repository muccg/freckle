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

	index=NULL;
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

// indexes the dots by y and then by x
void DotStore::CreateIndex()
{
	// if there is an old index destroy it
	if(index)
		DestroyIndex();

	// create a new index
	index=new LinkedListVal<IndexYNode *>;

	// go through every dot
	for(int i=0; i<numdots; i++)
		IndexDot(GetDot(i));
}

void DotStore::IndexDot( Dot *dot )
{
	assert(dot);

	// index this dot into the linkedlist index
	int x=dot->x;
	int y=dot->y;
	
	IndexXNode *xnode=new IndexXNode;
	memset(xnode,0,sizeof(IndexXNode)); 

	xnode->x=x;
	xnode->dot=dot;

	IndexYNode *ynode=NULL;

	// find this y in our linked list. loop through the y values
	LinkedListVal<IndexYNode*>::Iterator i(*index);
	for(; !i.Done() && (*i)->y<y; i++)
		;

	
	// if i is done, then we have scanned the whole list and every y is less than us. so insert a new y at the end of the list
	if(i.Done())
	{
		ynode=new IndexYNode;
		memset(ynode, 0, sizeof(IndexYNode));
		
		ynode->y=y;
		ynode->child=new LinkedListVal<IndexXNode*>;
		index->Append(ynode);
	}
	// if the iterators node is equal in y value to us, then we belong under this node
	else if((*i)->y==y)
	{
		ynode=*i;
	}
	// we need to be created and inserted after the present iterator point
	else
	{
		ynode=new IndexYNode;
		memset(ynode, 0, sizeof(IndexYNode));
		
		ynode->y=y;
		ynode->child=new LinkedListVal<IndexXNode*>;;
		index->AddAfter((LinkedListValElement<IndexYNode*>*)(*i), ynode);		// explicitly cast? WTF? Why?
	}

	assert(ynode);				//this should not be NULL by now. it should point to the ynode we belong under
	
	// now we work out where in this node we belong
	// if the node is empty then we can just add out xnode to the list
	if(!ynode->child->Size())
	{
		//add our xnode
		ynode->child->Append(xnode);
	}
	else
	{
		//find where we sit on the list
		LinkedListVal<IndexXNode *>::Iterator j(*(ynode->child));
		for(; !j.Done() && (*j)->x<x; j++)
			;

		if(j.Done())
		{
			//end of the list. append
			ynode->child->Append(xnode);
		}
		else
		{
			//mid list. insert
			ynode->child->AddAfter( (LinkedListValElement<IndexXNode*>*)(*j), xnode);
		}
	}
}

void DotStore::DestroyIndex()
{
	// nuke the index
	assert(index);	

	//for every y
	LinkedListVal<IndexYNode*>::Iterator i(*index);
	for(; !i.Done(); i++)
	{
		//for every x
		LinkedListVal<IndexXNode*>::Iterator j(*((*i)->child));
		for(; !j.Done(); j++)
		{
			// delete this xnode
			delete *j;

		}
		
		// delete this list
		(*i)->child->Empty();

		// delete this ynode
		delete *i;
	}
	
	// delete this list
	index->Empty();

	delete index;

	index=NULL;
}


//Dump the index to stdout. Mainly for DEBUG
void DotStore::DumpIndex()
{
	assert(index);
	printf("DumpIndex()\n");

	//for every y
	LinkedListVal<IndexYNode*>::Iterator i(*index);
	//printf("%d,%d,%d\n",(int)index->head,(int)index->tail,(int)index->Length());
	for(; !i.Done(); i++)
	{
		printf("y:%d\n",(*i)->y);
		//for every x
		LinkedListVal<IndexXNode*>::Iterator j(*((*i)->child));
		for(; !j.Done(); j++)
		{
			printf("\tx:%d => %d\n",(*j)->x,(int)(*j)->dot);
		}
	}
	printf("dumped\n");
}


// use the index to quickly find a dot
Dot *DotStore::GetIndexDot(int x,int y)
{
	assert(index);				//we must be indexed
	
	LinkedListVal<IndexYNode*>::Iterator i(*index);
	for(; !i.Done() && (*i)->y<y; i++)
		;
	
	// did we find the y
	if(i.Done() || (*i)->y!=y)
		return NULL;

	// search in this for the x.
	LinkedListVal<IndexXNode *>::Iterator j(*((*i)->child));
	for(; !j.Done() && (*j)->x<x; j++)
		;

	// did we find the x
	if(j.Done() || (*j)->x!=x)
		return NULL;

	// found it
	return (*j)->dot;
}

int DotStore::CountAreaMatches(int x1, int y1, int x2, int y2, int window)
{
	assert(index);				// we must be indexed
	assert(x2-x1 == y2-y1);			// we must be square. TODO: add support for rectangular area

	int count=0;

	// for each y value that is interesting
	LinkedListVal<IndexYNode*>::Iterator i(*index);
	for(; !i.Done() && (*i)->y<=y1-window; i++)
		;

	//if not done
	if(!i.Done())
		//now loop through each applicable y value
		for(;!i.Done() && (*i)->y<y2; i++)
		{
			// applicable y value.
			// go through each x.
			LinkedListVal<IndexXNode *>::Iterator j(*((*i)->child));
			for(; !j.Done() && (*j)->x<x1-window; j++)
				;

			if(!j.Done())
				for(;!j.Done() && (*j)->x<x2; j++)
				{
					// applicable x and y
					int x=(*j)->x;
					int y=(*i)->y;
					int length=(*j)->dot->length;
					int protrude=length;			//how much protrudes into this calculation square
					
					int initialcount=count;
// 					printf("%d,%d,%d,%d (%d)\n",x,y,length,protrude);

					// which zone are we in?
					// 1. inside the window
					if( x>=x1 && x<=x2 && y>=y1 && y<=y2)
					{
// 						printf("1\n");
						//we are inside the window
						// truncate protrude if we extend outside to the bottom or the right of the window
						if(x+protrude > x2)
							protrude=x2-x;
						if(y+protrude > y2)
							protrude=y2-x;
						
						//scan down and to the right to see when our next match point comes up or until our length is exhausted. add one for each point
						int xp=x;
						int yp=y;
						do
						{ 
							count++;
							protrude--;
						} while( (!GetIndexDot(++xp,++yp)) && protrude);
					}
					// 2. The parallelogram above the area
					else if( x>(y-y1+x1) && x<=(y-y1+x2) )
					{
// 						printf("2\n");
						//we are above the window
						if(length>y1-y)
						{
							//we extend into the window
							protrude=length-y+y1;
							
// 							printf("length:%d - y:%d + y1:%d = protrude:%d\n",length,y,y1,protrude);
	
							// check if we extend out of the window to the right too
							int sigma=length-x2+x;
							if(sigma>0)
								protrude-=sigma;

							// TODO: Check if we extend out of the bottom (non square window)
// 							printf("sigma:%d protrude:%d\n",sigma,protrude);
	
							int xp=x;
							int yp=y;
							do
							{
								if(yp>=y1)
								{
									count++;
									protrude--;
								}
							} while((!GetIndexDot(++xp,++yp)) && protrude);
						}
					}
					// 3. The parallelogram to the left of the area
					else if( y>=(x-x1+y1) && y<=(x-x1+y2) )
					{
// 						printf("3\n");
						//we are to the left of the window
						if(length>x1-x)
						{
							protrude=length-x+x1;
							
							// check if we extend out of the window to the bottom too
							int sigma=length-y2+y;
							if(sigma>0)
								protrude-=sigma;
	
							// TODO: Check if we extend out of the window to the right, too (non square window)
	
							int xp=x;
							int yp=y;
							do
							{
								if(xp>=x1)
								{
									count++;
									protrude--;
								}
							} while((!GetIndexDot(++xp,++yp)) && protrude);
						}
					}
					// else we are not included.

// 					printf("count=%d\n",count-initialcount);
				}
		}

// 	printf("count=%d\n",count);
	return count;
}


