#include "QuadTreeNode.h"
#include <memory.h>

QuadTreeNode::QuadTreeNode()
{
	// zero everything
	memset(this, 0, sizeof(QuadTreeNode));
	type=QTNODE;
}

QuadTreeNode::QuadTreeNode(Type xp, Type yp, Type xp1, Type yp1, Type xp2, Type yp2)
{
	store.child[NW]=NULL;
	store.child[NE]=NULL;
	store.child[SW]=NULL;
	store.child[SE]=NULL;

	// our division point	
	x=xp;
	y=yp;

	// our inclusinve extents
	x1=xp1;
	y1=yp1;
	x2=yp2;
	y2=xp2;

	type=QTNODE;
}

// this takes the extents of the node area and calculates the 1/2 way point and does everything needed
// THIS MAKES IT A LEAF NODE
QuadTreeNode::QuadTreeNode(Type xp1, Type yp1, Type xp2, Type yp2)
{
	// zero everything
	memset(this, 0, sizeof(QuadTreeNode));

	x1=xp1;
	x2=xp2;
	y1=yp1;
	y2=yp2;
	type=QTLEAF;
}

QuadTreeNode::~QuadTreeNode()
{
	if(type==QTNODE)
	{
		for(int direction=0; direction<4; direction++)
			if(store.child[direction])
				delete store.child[direction];
	}
}

void QuadTreeNode::AddDot(Dot *dt)
{
	if(isNode())
	{
		if(dt->x < x)
		{
			if(dt->y < y)
			{	
				if(!store.child[NW])
					store.child[NW]=new QuadTreeNode(x1,y1,x,y);
				store.child[NW]->AddDot(dt);
			}
			else
			{
				if(!store.child[SW])
					store.child[SW]=new QuadTreeNode(x1,y,x,y2);
				store.child[SW]->AddDot(dt);
			}
		}
		else
		{
			if(dt->y < y)
			{
				if(!store.child[NE])
					store.child[NE]=new QuadTreeNode(x,y1,x2,y);
				store.child[NE]->AddDot(dt);
			}
			else
			{
				if(!store.child[SE])
					store.child[SE]=new QuadTreeNode(x,y,x2,y2);
				store.child[SE]->AddDot(dt);
			}
		}
	}
	else
	{
		// this is a leaf. see if we can store it here
		for(int i=0; i<NUMDOTS; i++)
			if(!store.dot[i])
			{
				store.dot[i]=dt;
				return;
			}
		
		// cant fit in this leaf. Lets turn the leaf to a node and then try re-adding it
		LeafToNode();
		AddDot(dt);	
	}

}

// Turn this leaf into a node and reclassify the contents
void QuadTreeNode::LeafToNode()
{
	assert(isLeaf());

	// save the dots
	Dot *dotstore[NUMDOTS];
	for(int i=0; i<NUMDOTS; i++)
		dotstore[i]=store.dot[i];
	
	// now change to a node and subdivide
	MakeNode((x2-x1)/2+x1, (y2-y1)/2+y1);

	for(int i=0; i<NUMDOTS; i++)
		AddDot(dotstore[i]);
	
}

// macro to determine if rectangle r1 and r2 intersect
// r1 is defined as (r1l,r1t)-(r1r,r1b)
// r2 is defined as (r2l,r2t)-(r2r,r2b)
#define RectIntersect(r1l,r1t,r1r,r1b,r2l,r2t,r2r,r2b)	\
	( ! ( r2l>r1r || r2r<r1l || r2t>r1b || r2b<r1t) )

// Issues a spatial query on the given rectangle. Assembled a linked list of the results and returns it
void QuadTreeNode::SpatialQueryRecurse(LinkedListVal<Dot *> *list, Type xp1, Type yp1, Type xp2, Type yp2)
{
	if(x1>xp2 || x2<xp1 || y1>yp2 || y2<yp1)
		return;

	if(isLeaf())
	{
		// check each dot and if its included, add it to the list
		for(int i=0; i<NUMDOTS; i++)
			if(store.dot[i])
				if(store.dot[i]->x >= xp1 && store.dot[i]->x <= xp2 && store.dot[i]->y >= yp1 && store.dot[i]->y <= yp2)
					list->Append(store.dot[i]);
		return;
	}
	else if(isNode())
	{
		if(store.child[NW] && RectIntersect(xp1,yp1,xp2,yp2,x1,y1,x,y))
			store.child[NW]->SpatialQueryRecurse(list, xp1, yp1, xp2, yp2);
		
		if(store.child[NE] && RectIntersect(xp1,yp1,xp2,yp2,x,y1,x2,y))
			store.child[NE]->SpatialQueryRecurse(list, xp1, yp1, xp2, yp2);

		if(store.child[SW] && RectIntersect(xp1,yp1,xp2,yp2,x1,y,x,y2))
			store.child[SW]->SpatialQueryRecurse(list, xp1, yp1, xp2, yp2);

		if(store.child[SE] && RectIntersect(xp1,yp1,xp2,yp2,x,y,x2,y2))
			store.child[SE]->SpatialQueryRecurse(list, xp1, yp1, xp2, yp2);
		
	}

}

// deletes a Dot from the store. Recursive function
void QuadTreeNode::DeleteDot(Dot *dot)
{
	assert(dot);

	int xp=dot->x;
	int yp=dot->y;

	if(isLeaf())
	{
		// check each dot in the leaf to see if we are one.
		for(int i=0; i<NUMDOTS; i++)
			if(store.dot[i])
				if(store.dot[i]->x ==xp && store.dot[i]->y == yp)
				{
					// this is the dot. lets delete it
					// we move the dots above down and set the end one to NULL
					for(int j=i; j<NUMDOTS-1;j++)
						store.dot[j]=store.dot[j+1];
					store.dot[NUMDOTS-1]=NULL;

					return;
				}
		
		// if we get to here then the Dot was not found in the tree, or the dot has moved since it was added
		// this is very bad and we should die
		Dump();
		assert(0);
	}
	else if(isNode())
	{
		if(store.child[NW] && xp<x && yp<y)
			store.child[NW]->DeleteDot(dot);
		
		if(store.child[NE] && xp>=x && yp<y)
			store.child[NE]->DeleteDot(dot);
		
		if(store.child[SW] && xp<x && yp>=y)
			store.child[SW]->DeleteDot(dot);
		
		if(store.child[SE] && xp>=x && yp>=y)
			store.child[SE]->DeleteDot(dot);
	}
}


