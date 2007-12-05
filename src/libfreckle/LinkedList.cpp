#include "LinkedList.h"
#include "assert.h"

template <class T>
LinkedList<T>::LinkedList()
{
	head=NULL;
	tail=NULL;
	num=0;
}

template <class T>
LinkedList<T>::~LinkedList()
{
	Empty();
}

template <class T>
void LinkedList<T>::AddItem(T *item)
{
	LLNode<T> *node=new LLNode<T>();
	node->SetData(item);

	if(!head || !tail || !num)
	{
		//first item. assert that we are consistent
		assert(!head);
		assert(!tail);
		assert(!num);

		head=tail=node;
		num=1;
	}
	else
	{
		// append to the end of the list
		tail->SetNext(node);
		node->SetPrev(tail);
		tail=node;
		num++;
	}	
}

template <class T>
void LinkedList<T>::DelItem(T *item)
{
	//find item
	LLNode<T> *node=head;
	for( ; node && node->GetData()!=item; node=node->GetNext() )
		;

	assert(node);				// not found!

	// close gap. 3 possible conditions. We are the head. We are the tail. We are both the head and the tail. Or we are somehwere in between
	if(head==node)
	{
		if(tail==node)
		{
			//only item on list
			head=tail=NULL;
		}
		else
		{
			//head. pull from start
			head=node->GetNext();
			head->SetPrev(NULL);
		}
	}
	else if(tail==node)
	{
		//tail. pop off end
		tail=node->GetPrev();
		tail->SetNext(NULL);
	}
	else
	{
		//middle. close gap
		node->GetPrev()->SetNext( node->GetNext() );
		node->GetNext()->SetPrev( node->GetPrev() );
	}
	num--;
	delete node;
}

template <class T>
int LinkedList<T>::IndexOf(T *item)
{
	int i=0;

	//find item
	for(LLNode<T> *node=head; node && node->GetData()!=item; node=node->GetNext() )
		i++;

	assert(i<num);

	return i;
}

template <class T>
T *LinkedList<T>::GetItem(int i)
{
	assert(i>=0);
	assert(i<num);

	LLNode<T> *node=head;
	for(int j=0; i<i; j++)
	{
		assert(node);
		node=node->GetNext();
	}

	return node->GetData();
}

template <class T>
int LinkedList<T>::GetNum()
{
	return num;
}

template <class T>
void LinkedList<T>::Empty()
{
	if(!head)
		return;				// already empty

	LLNode<T> *nextnode=NULL;
	for(LLNode<T> *node=head; node; node=nextnode )
	{
		nextnode=node->GetNext();
		delete node;
	}
	num=0;
	head=NULL;
	tail=NULL;
}






