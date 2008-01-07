#ifndef _QUADTREE_H_
#define _QUADTREE_H_

#include "QuadTreeNode.h"
#include "LinkedListVal.h"
#include <stdio.h>

class QuadTree
{
private:
	QuadTreeNode	*root;

	Type		minx,miny,maxx,maxy;			// The bounds of our world
	Type		width, height;

public:
	QuadTree(Type x1, Type y1, Type x2, Type y2);
	~QuadTree();

	//! \brief add a dots location to the index 
	void AddDot(Dot *dot);

	//! \brief delete a dot from the index. Dot must have the same position as where it sits in the tree
	void DelDot(Dot *dot);				// dot MUST not have moved from its insertion position

	LinkedListVal<Dot *> *SpatialQuery(Type x1, Type y1, Type x2, Type y2);

	// dump the tree to stdout
	inline void Dump()
	{
		printf("QuadTree(%d)\n=======================\n",(int)this);
		printf("minx:%d\tminy:%d\n",minx,miny);
		printf("maxx:%d\tmaxy:%d\n",maxx,maxy);
		printf("width:%d\theight:%d\n",width,height);
		printf("\n");
		if(root)
			root->Dump();
	}

};










#endif

