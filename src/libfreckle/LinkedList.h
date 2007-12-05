#ifndef _LINKEDLIST_H_
#define _LINKEDLIST_H_

#include <iostream>
#define LIST_IS_EMPTY 1

template<class T>
class LLNode
{
private:
	// pointers up and down the list
	LLNode<T>	*next;
	LLNode<T> *prev;

	// our item
	T	*data;

public:
	LLNode();
	~LLNode();

	inline void SetNext(LLNode<T> *n)
	{
		next=n;
	}
	
	inline void SetPrev(LLNode<T> *p)
	{
		prev=p;
	}
	
	inline void SetData(T *item)
	{
		data=item;
	}
	
	inline LLNode<T> *GetNext()
	{
		return next;
	}
	
	inline LLNode<T> *GetPrev()
	{
		return prev;
	}

	inline T *GetData()
	{
		return data;
	}
};

template<class T>
class LinkedList
{
private:
	// head and tail
	LLNode<T> *head;
	LLNode<T> *tail;
	
	int num;

public:
	LinkedList();
	~LinkedList();

	// add and remove items
	void AddItem(T *item);
	void DelItem(T *item);

	// find and items index
	int IndexOf(T *item);
	
	// get item at index
	T *GetItem(int i);

	// get number of items
	int GetNum();

	// empty the list
	void Empty();
};


#endif
