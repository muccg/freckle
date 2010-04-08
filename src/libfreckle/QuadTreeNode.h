#ifndef _QUADTREENODE_H_
#define _QUADTREENODE_H_

#include "Dot.h"
#include "LinkedListVal.h"
#include <assert.h>
#include <stdio.h>
#include <types.h>

// the type of units for the points
typedef int Type;

#define NW	0
#define NE	1
#define SW	2
#define SE	3

#define QTLEAF	0
#define QTNODE	1

#define NUMDOTS 16


class QuadTreeNode
{
private:
	union
	{
		QuadTreeNode	*child[4];				// if we are a node
		Dot		*dot[NUMDOTS];				// if we are a leaf, and could contain up to 4 dots
	} store;

	Type		x,y;						// the centre point of the split (for non uniform splitting)
	int		type;						// leaf or node

	// the extents of this node (<=, >=)
	Type 		x1,y1,x2,y2;

public:
	QuadTreeNode();
	QuadTreeNode(Type xp, Type yp, Type xp1, Type yp1, Type xp2, Type yp2);
	QuadTreeNode(Type x1, Type y1, Type x2, Type y2);
	~QuadTreeNode();

	//! \brief recursively add a dot
	void AddDot(Dot *dot);

	//! \brief deletes a Dot from the store. Recursive function
	void DeleteDot(Dot *dot);

	// turn our existing leaf into a node
	void LeafToNode();

	inline bool isNode() const
	{
		return type==QTNODE;
	}

	inline bool isLeaf() const
	{
		return type==QTLEAF;
	}

	inline void MakeNode(Type xp, Type yp)
	{
		type=QTNODE;

		x=xp;
		y=yp;

		// no nodes yet. they are created on demand
		store.child[NW]=NULL;
		store.child[NE]=NULL;
		store.child[SW]=NULL;
		store.child[SE]=NULL;
	}

	inline void MakeLeaf()
	{
		type=QTLEAF;

		x=y=0;

		// TODO recursilvely delete any child nodes
	}

	inline Dot *GetDot(int index) const
	{
		assert(index>=0 && index<NUMDOTS);
		assert(isLeaf());
		return store.dot[index];
	}

	inline QuadTreeNode *GetChild(int direction) const
	{
		assert(direction>=0 && direction<4);
		assert(isNode());
		return store.child[direction];
	}

	inline void SetDot(int index, Dot *dot)
	{
		assert(index>=0 && index<NUMDOTS);
		assert(isLeaf());
		store.dot[index]=dot;
	}

	inline void SetChild(int direction, QuadTreeNode *child)
	{
		assert(direction>=0 && direction<4);
		assert(isNode());
		store.child[direction]=child;
	}

	inline Type GetX() const
	{
		return x;
	}

	inline Type GetY() const
	{
		return y;
	}

	inline void SetX(Type xp)
	{
		x=xp;
	}

	inline void SetY(Type yp)
	{
		y=yp;
	}

	inline void SetXY(Type xp, Type yp)
	{
		x=xp;
		y=yp;
	}
	
#ifdef DEBUG
	inline void Dump()
	{
		printf("QuadTreeNode(%d)\n============================\n",(__U32_TYPE)this);
		printf("type:%s\n",type==QTNODE?"QTNODE":"QTLEAF");
		printf("x1:%d\ty1:%d\tx2:%d\ty2:%d\n",x1,y1,x2,y2);
		if(type==QTNODE)
		{
			printf("x:%d\ty:%d\n",x,y);
			printf("children: %d %d %d %d\n\n",(__U32_TYPE)store.child[0],(__U32_TYPE)store.child[1],(__U32_TYPE)store.child[2],(__U32_TYPE)store.child[3]);
			for(int i=0; i<4; i++)
				if(store.child[i])
					store.child[i]->Dump();
		}
		else
		{
			printf("dots:\n");
			for(int i=0; i<NUMDOTS; i++)
				if(store.dot[i])
					printf("dot %d: %d,%d,%d\n",(int)store.dot[i],store.dot[i]->x,store.dot[i]->y,store.dot[i]->length);
			printf("\n");
		}
	}
#endif

	void SpatialQueryRecurse(LinkedListVal<Dot *> *list, Type xp1, Type yp1, Type xp2, Type yp2);

};













#endif

