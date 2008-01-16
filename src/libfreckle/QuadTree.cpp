#include "QuadTree.h"
#include "LinkedListVal.h"

#include <assert.h>
#include <stdio.h>

QuadTree::QuadTree(Type x1, Type y1, Type x2, Type y2)
{
	root=NULL;

	printf("x1:%d, y1:%d, x2:%d, y2:%d\n", x1, y1, x2, y2);

	assert(x1<x2);
	assert(y1<y2);

	minx=x1;
	miny=y1;
	maxx=x2;
	maxy=y2;

	width=x2-x1;
	height=y2-y1;
}

QuadTree::~QuadTree()
{
	if(root)
		delete root;
}


void QuadTree::AddDot(Dot *dot)
{
	assert(dot);
	if(!root)
	{
		//first quad node
		root=new QuadTreeNode(minx,miny,maxx,maxy);
	}

	// Add to our root
	root->AddDot(dot);
}

void QuadTree::DelDot(Dot *dot)
{
	assert(dot);
	assert(root);

	//delete the dot from the quadtree
	root->DeleteDot(dot);
	
}

// Issues a spatial query on the given rectangle. Assembled a linked list of the results and returns it
LinkedListVal<Dot *> *QuadTree::SpatialQuery(Type xp1, Type yp1, Type xp2, Type yp2)
{
	assert(xp1<=xp2);
	assert(yp1<=yp2);
	
	LinkedListVal<Dot *> *list=new LinkedListVal<Dot *>;

	// walk the tree looking for whats potentially included and whats not
	if(root)							// if root is NULL, then our index is empty		
		root->SpatialQueryRecurse(list,xp1,yp1,xp2,yp2);

	return list;
}
