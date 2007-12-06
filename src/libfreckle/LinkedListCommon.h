// LinkedList
// ==========
// TODO: rewrite
//
// $Id: LinkedListCommon.hh,v 1.8.2.1 2003/07/29 12:27:37 crispin Exp $
// $Revision: 1.8.2.1 $

#if defined( LINKEDLIST_VALUE )
/**
 * \brief Double list element for LinkedListVal.
 *
 * LinkedListValElement stores the type T of LinkedListVal by value. 
 */
template <class T> struct LinkedListValElement
{
	/**
	 * \brief Construct a LinkedListValElement with a given value.
	 *
	 * The only constructor available initializes the value element. This
	 * enforces that LinkedListVal elements are never created without having their
	 * value intialzed by the user. T's copy constructor is used to copy the
	 * value in.
	 */
	LinkedListValElement( const T &val ) : value(val) { }

	/**
	 * \brief Value stored by the list element.
	 *
	 * Value is always copied into new list elements using the copy
	 * constructor.
	 */
	T value;

	/**
	 * \brief List previous pointer.
	 *
	 * Points to the previous item in the list. If this is the First item in
	 * the list, then prev is NULL. If this element is not in a list then
	 * prev is undefined.
	 */
	LinkedListValElement<T> *prev;

	/**
	 * \brief List next pointer.
	 *
	 * Points to the next item in the list. If this is the list item in the
	 * list, then next is NULL. If this element is not in a list then next is
	 * undefined.
	 */
	LinkedListValElement<T> *next;
};
#else

#ifndef DOUBLE_LIST_EL
#define DOUBLE_LIST_EL
/**
 * \brief Double list element properties.
 *
 * This class can be inherited to make a class suitable to be a double list
 * element. It simply provides the next and previous pointers. An alternative
 * is to put the next and previous pointers in the class directly.
 */
template <class Element> struct LinkedListElement
{
	/**
	 * \brief List previous pointer.
	 *
	 * Points to the previous item in the list. If this is the First item in
	 * the list, then prev is NULL. If this element is not in a list then
	 * prev is undefined.
	 */
	Element *prev;

	/**
	 * \brief List next pointer.
	 *
	 * Points to the next item in the list. If this is the list item in the
	 * list, then next is NULL. If this element is not in a list then next is
	 * undefined.
	 */
	Element *next;
};
#endif /* DOUBLE_LIST_EL */

#endif

/* Doubly Linked List */
template <LINKEDLIST_MULTIELEMENT_TEMPDEF> class LinkedList
{
public:
	/** \brief Initialize an Empty list. */
	LinkedList() : head(0), tail(0), listLen(0) {}

	LinkedList(const LinkedList &other);

#ifdef LINKEDLIST_VALUE
	/**
	 * \brief Clear the double list contents.
	 *
	 * All elements are deleted.
	 */
	~LinkedList() { Empty(); }
#else
	/**
	 * \brief Abandon all elements in the list. 
	 *
	 * List elements are not deleted.
	 */
	~LinkedList() {}
#endif


	/* Shallow and deep copy. */
	void shallowCopy(const LinkedList &other);

	/* Assignment operator. */
	LinkedList &operator=(const LinkedList &other);

#ifdef LINKEDLIST_VALUE
	/**
	 * \brief Make a new element and Prepend it to the front of the list.
	 *
	 * The item is copied into the new element using the copy constructor.
	 * Equivalent to list.AddBefore(list.head, item).
	 */
	void Prepend(const T &item);

	/**
	 * \brief Make a new element and Append it to the end of the list.
	 *
	 * The item is copied into the new element using the copy constructor.
	 * Equivalent to list.AddAfter(list.tail, item).
	 */
	void Append(const T &item);

	/**
	 * \brief Make a new element and insert it immediately after an element in
	 * the list.
	 *
	 * The item is copied into the new element using the copy constructor. If
	 * prev_el is NULL then the new element is Prepended to the front of the
	 * list. If prev_el is not already in the list then undefined behaviour
	 * results.  Equivalent to list.AddAfter(prev_el, new LinkedListValElement(item)).
	 */
	void AddAfter(Element *prev_el, const T &item);

	/**
	 * \brief Make a new element and insert it immediately before an element
	 * in the list. 
	 *
	 * The item is copied into the new element using the copy construcotor. If
	 * next_el is NULL then the new element is Appended to the end of the
	 * list.  If next_el is not already in the list then undefined behaviour
	 * results.  Equivalent to list.AddBefore(next_el, new LinkedListValElement(item)).
	 */
	void AddBefore(Element *next_el, const T &item);
#endif

	/**
	 * \brief Prepend a single element to the front of the list.
	 *
	 * If new_el is already an element of some list, then undefined behaviour
	 * results. Equivalent to list.AddBefore(list.head, new_el).
	 */
	void Prepend(Element *new_el) { AddBefore(head, new_el); }

	/**
	 * \brief Append a single element to the end of the list.
	 *
	 * If new_el is alreay an element of some list, then undefined behaviour
	 * results.  Equivalent to list.AddAfter(list.tail, new_el).
	 */
	void Append(Element *new_el)  { AddAfter(tail, new_el); }

	/**
	 * \brief Prepend an entire list to the beginning of this list.
	 *
	 * All items are moved, not copied. Afterwards, the other list is emtpy.
	 * All items are Prepended at once, so this is an O(1) operation.
	 * Equivalent to list.AddBefore(list.head, dl).
	 */
	void Prepend(LinkedList &dl)       { AddBefore(head, dl); }

	/**
	 * \brief Append an entire list to the end of the list.
	 *
	 * All items are moved, not copied. Afterwards, the other list is Empty.
	 * All items are appened at once, so this is an O(1) operation.
	 * Equivalent to list.AddAfter(list.tail, dl).
	 */
	void Append(LinkedList &dl)        { AddAfter(tail, dl); }

	void AddAfter(Element *prev_el, Element *new_el);
	void AddBefore(Element *next_el, Element *new_el);

	void AddAfter(Element *prev_el, LinkedList &dl);
	void AddBefore(Element *next_el, LinkedList &dl);

	/**
	 * \brief Detach the head of the list
	 *
	 * The element Detached is not deleted. If there is no head of the list
	 * (the list is Empty) then undefined behaviour results.  Equivalent to
	 * list.Detach(list.head).
	 *
	 * \returns The element Detached.
	 */
	Element *DetachFirst()        { return Detach(head); }

	/**
	 * \brief Detach the tail of the list
	 *
	 * The element Detached is not deleted. If there is no tail of the list
	 * (the list is Empty) then undefined behaviour results.  Equivalent to
	 * list.Detach(list.tail).
	 *
	 * \returns The element Detached.
	 */
	Element *DetachLast()         { return Detach(tail); }

 	/* Detaches an element from the list. Does not free any memory. */
	Element *Detach(Element *el);

	//! pop the last element off the list
	//! deletes the container as it does
	LINKEDLIST_MULTIELEMENT_TEMPUSE Pop()				
	{ 
		Element *result=DetachLast(); 
		LINKEDLIST_MULTIELEMENT_TEMPUSE answer=result->value; 
		delete result; 
		return answer;
	}
	
	//! pop the first element off the list
	//! deletes the container as it does
	LINKEDLIST_MULTIELEMENT_TEMPUSE PrePop()				
	{ 
		Element *result=DetachFirst(); 
		LINKEDLIST_MULTIELEMENT_TEMPUSE answer=result->value; 
		delete result; 
		return answer;
	}
	
	/**
	 * \brief Detach and delete the First element in the list.
	 *
	 * If there is no First element (the list is Empty) then undefined
	 * behaviour results.  Equivalent to delete list.Detach(list.head);
	 */
	void RemoveFirst()         { delete Detach( head ); }

	/**
	 * \brief Detach and delete the Last element in the list.
	 *
	 * If there is no Last element (the list is emtpy) then undefined
	 * behaviour results.  Equivalent to delete list.Detach(list.tail);
	 */
	void RemoveLast()          { delete Detach( tail ); }

	/**
	 * \brief Detach and delete an element from the list.
	 *
	 * If the element is not in the list, then undefined behaviour results.
	 * Equivalent to delete list.Detach(el);
	 */
	void Remove(Element *el)   { delete Detach( el ); }

	//
	// Find(LINKEDLIST_MULTIELEMENT_TEMPUSE val)
	//
	// return the element associated with the value
	//
	Element *Find(LINKEDLIST_MULTIELEMENT_TEMPUSE val)
	{
		for(Element *point=head; point; point=point->next)
			if(point->value==val)
				return point;				
		
		return (Element *)0;
	}
	
#ifdef LINKEDLIST_VALUE
	T &operator[](int a)
	{
		int j=0;
				
		Element *point=head;
		for(; j<a; j++)
			point=point->next;
		
		return point->value;
	}
#endif
	
	
	//! Remove an element from the list by passed value
	void Remove(LINKEDLIST_MULTIELEMENT_TEMPUSE val)
	{
		Remove(Find(val));
	}
	
	void Empty();
	void Abandon();

	/**
	 * \brief Return the length.
	 *
	 * \returns The number of elements in the list. 
	 */
	int Length() const { return listLen; }

	/** \brief Head and tail of the linked list. */
	Element *head, *tail;

	/** \brief The number of element in the list. */
	int listLen;

	/* Convenience access. */
	int Size() const           { return listLen; }

	/* Forward this so a ref can be used. */
	struct Iterator;

	/* Class for setting the iterator. */
	struct IteratorFirst { IteratorFirst( const LinkedList &l ) : l(l) { } const LinkedList &l; };
	struct IteratorLast { IteratorLast( const LinkedList &l ) : l(l) { } const LinkedList &l; };
	struct IteratorNext { IteratorNext( const Iterator &i ) : i(i) { } const Iterator &i; };
	struct IteratorPrev { IteratorPrev( const Iterator &i ) : i(i) { } const Iterator &i; };

	/* Double list iterator. */
	struct Iterator
	{
		/* Default construct. */
		Iterator() : ptr(0) { }

		/* Construct from a double list. */
		Iterator( const LinkedList &dl )      : ptr(dl.head) { }
		Iterator( const IteratorFirst &dlf ) : ptr(dlf.l.head) { }
		Iterator( const IteratorLast &dll )  : ptr(dll.l.tail) { }
		Iterator( const IteratorNext &dln )  : ptr(dln.i.ptr->BASEREFERENCE(next)) { }
		Iterator( const IteratorPrev &dlp )  : ptr(dlp.i.ptr->BASEREFERENCE(prev)) { }

		/* Assign from a double list. */
		Iterator &operator=( const LinkedList &dl )     { ptr = dl.head; return *this; }
		Iterator &operator=( const IteratorFirst &af ) { ptr = af.l.head; return *this; }
		Iterator &operator=( const IteratorLast &al )  { ptr = al.l.tail; return *this; }
		Iterator &operator=( const IteratorNext &an )  { ptr = an.i.ptr->BASEREFERENCE(next); return *this; }
		Iterator &operator=( const IteratorPrev &ap )  { ptr = ap.i.ptr->BASEREFERENCE(prev); return *this; }

		/* At the end, beginning? */
		bool More() const    { return ptr != 0; }
		bool Done() const    { return ptr == 0; }
		bool RevMore() const { return ptr != 0; }
		bool RevDone() const { return ptr == 0; }

		/* At the First, Last element. */
		bool First() const { return ptr && ptr->BASEREFERENCE(prev) == 0; }
		bool Last() const  { return ptr && ptr->BASEREFERENCE(next) == 0; }

#ifdef LINKEDLIST_VALUE
		/* Cast, dereference, arrow ops. */
		operator T*() const         { return &ptr->value; }
		T &operator *() const       { return ptr->value; }
		T *operator->() const       { return &ptr->value; }
#else
		/* Cast, dereference, arrow ops. */
		operator Element*() const   { return ptr; }
		Element &operator *() const { return *ptr; }
		Element *operator->() const { return ptr; }
#endif

		/* Increment. */
		inline Element *operator++()      { return ptr = ptr->BASEREFERENCE(next); }
		inline Element *increment()       { return ptr = ptr->BASEREFERENCE(next); }
		inline Element *operator++(int);

		/* Decrement. */
		inline Element *operator--()      { return ptr = ptr->BASEREFERENCE(prev); }
		inline Element *decrement()       { return ptr = ptr->BASEREFERENCE(prev); }
		inline Element *operator--(int);

		/* Return the next, prev. Does not modify. */
		inline IteratorNext next() const { return IteratorNext(*this); }
		inline IteratorPrev prev() const { return IteratorPrev(*this); }

		/* The iterator is simply a pointer. */
		Element *ptr;
	};

	/* Return classes for setting First, Last. */
	IteratorFirst First()  { return IteratorFirst(*this); }
	IteratorLast Last()    { return IteratorLast(*this); }
};

/** 
 * \brief Perform a deep copy of the list.
 * 
 * The elements of the other list are duplicated and put into this list.
 * Elements are copied using the copy constructor.
 */
template <LINKEDLIST_MULTIELEMENT_TEMPDEF> LinkedList<LINKEDLIST_MULTIELEMENT_TEMPUSE>::
		LinkedList(const LinkedList<LINKEDLIST_MULTIELEMENT_TEMPUSE> &other) :
			head(0), tail(0), listLen(0)
{
	Element *el = other.head;
	while( el != 0 ) {
		Append( new Element(*el) );
		el = el->BASEREFERENCE(next);
	}
}

/**
 * \brief Shallow copy another list into this list.
 *
 * The elements of the other list are copied in by reference. The two lists
 * will share the same list elements. If this list contains any elements
 * before the copy, then they are Abandoned. If either of the lists are
 * modified after a shallow copy, the other list may become corrupted.
 */
template <LINKEDLIST_MULTIELEMENT_TEMPDEF> void LinkedList<LINKEDLIST_MULTIELEMENT_TEMPUSE>::
		shallowCopy(const LinkedList &other)
{
	head = other.head;
	listLen = other.listLen;
}

#ifdef LINKEDLIST_VALUE

/**
 * \brief Assign another list into this list using a deep copy.
 *
 * The elements of the other list are duplicated and put into this list.  Each
 * list item is created using the copy constructor. If this list contains any
 * elements before the copy, they are deleted First.
 *
 * \returns A refence to this.
 */
template <LINKEDLIST_MULTIELEMENT_TEMPDEF> LinkedList<LINKEDLIST_MULTIELEMENT_TEMPUSE> &LinkedList<LINKEDLIST_MULTIELEMENT_TEMPUSE>::
		operator=(const LinkedList &other)
{
	/* Loose the old list. The value list assumes items were allocated on the
	 * stack by ourselves. It assumes ownership. */
	Empty();

	Element *el = other.head;
	while( el != 0 ) {
		Append( new Element(*el) );
		el = el->BASEREFERENCE(next);
	}
	return *this;
}

/* Prepend a new item. Inlining this bloats the caller with new overhead. */
template <LINKEDLIST_MULTIELEMENT_TEMPDEF> void LinkedList<LINKEDLIST_MULTIELEMENT_TEMPUSE>::
		Prepend(const T &item)
{
	AddBefore(head, new Element(item)); 
}

/* Append a new item. Inlining this bloats the caller with the new overhead. */
template <LINKEDLIST_MULTIELEMENT_TEMPDEF> void LinkedList<LINKEDLIST_MULTIELEMENT_TEMPUSE>::
		Append(const T &item)
{
	AddAfter(tail, new Element(item));
}

/* Add a new item after a prev element. Inlining this bloats the caller with
 * the new overhead. */
template <LINKEDLIST_MULTIELEMENT_TEMPDEF> void LinkedList<LINKEDLIST_MULTIELEMENT_TEMPUSE>::
		AddAfter(Element *prev_el, const T &item)
{
	AddAfter(prev_el, new Element(item));
}

/* Add a new item before a next element. Inlining this bloats the caller with
 * the new overhead. */
template <LINKEDLIST_MULTIELEMENT_TEMPDEF> void LinkedList<LINKEDLIST_MULTIELEMENT_TEMPUSE>::
		AddBefore(Element *next_el, const T &item)
{
	AddBefore(next_el, new Element(item));
}

#else

/**
 * \brief Assign another list into this list using a deep copy.
 *
 * The elements of the other list are duplicated and put into this list.  Each
 * list item is created using the copy constructor. If this list contains any
 * elements before the copy, they are Abandoned.
 *
 * \returns A refence to this.
 */
template <LINKEDLIST_MULTIELEMENT_TEMPDEF> LinkedList<LINKEDLIST_MULTIELEMENT_TEMPUSE> &LinkedList<LINKEDLIST_MULTIELEMENT_TEMPUSE>::
		operator=(const LinkedList &other)
{
	/* Non value lists do not assume ownership of the data. The old data is
	 * simply lost. */
	head = 0;
	listLen = 0;

	Element *el = other.head;
	while( el != 0 ) {
		Append( new Element(*el) );
		el = el->BASEREFERENCE(next);
	}
	return *this;
}

#endif

/*
 * The larger iterator operators.
 */

/* Postfix ++ */
template <LINKEDLIST_MULTIELEMENT_TEMPDEF> Element *LinkedList<LINKEDLIST_MULTIELEMENT_TEMPUSE>::Iterator::
		operator++(int)       
{
	Element *rtn = ptr; 
	ptr = ptr->BASEREFERENCE(next);
	return rtn;
}

/* Postfix -- */
template <LINKEDLIST_MULTIELEMENT_TEMPDEF> Element *LinkedList<LINKEDLIST_MULTIELEMENT_TEMPUSE>::Iterator::
		operator--(int)       
{
	Element *rtn = ptr;
	ptr = ptr->BASEREFERENCE(prev);
	return rtn;
}

/**
 * \brief Insert an element immediately after an element in the list.
 *
 * If prev_el is NULL then new_el is Prepended to the front of the list. If
 * prev_el is not in the list or if new_el is already in a list, then
 * undefined behaviour results.
 */
template <LINKEDLIST_MULTIELEMENT_TEMPDEF> void LinkedList<LINKEDLIST_MULTIELEMENT_TEMPUSE>::
		AddAfter(Element *prev_el, Element *new_el)
{
	/* Set the previous pointer of new_el to prev_el. We do
	 * this regardless of the state of the list. */
	new_el->BASEREFERENCE(prev) = prev_el; 

	/* Set forward pointers. */
	if (prev_el == 0) {
		/* There was no prev_el, we are inserting at the head. */
		new_el->BASEREFERENCE(next) = head;
		head = new_el;
	} 
	else {
		/* There was a prev_el, we can access previous next. */
		new_el->BASEREFERENCE(next) = prev_el->BASEREFERENCE(next);
		prev_el->BASEREFERENCE(next) = new_el;
	} 

	/* Set reverse pointers. */
	if (new_el->BASEREFERENCE(next) == 0) {
		/* There is no next element. Set the tail pointer. */
		tail = new_el;
	}
	else {
		/* There is a next element. Set it's prev pointer. */
		new_el->BASEREFERENCE(next)->BASEREFERENCE(prev) = new_el;
	}

	/* Update list length. */
	listLen++;
}

/**
 * \brief Insert an element immediatly before an element in the list.
 *
 * If next_el is NULL then new_el is Appended to the end of the list. If
 * next_el is not in the list or if new_el is already in a list, then
 * undefined behaviour results.
 */
template <LINKEDLIST_MULTIELEMENT_TEMPDEF> void LinkedList<LINKEDLIST_MULTIELEMENT_TEMPUSE>::
		AddBefore(Element *next_el, Element *new_el)
{
	/* Set the next pointer of the new element to next_el. We do
	 * this regardless of the state of the list. */
	new_el->BASEREFERENCE(next) = next_el; 

	/* Set reverse pointers. */
	if (next_el == 0) {
		/* There is no next elememnt. We are inserting at the tail. */
		new_el->BASEREFERENCE(prev) = tail;
		tail = new_el;
	} 
	else {
		/* There is a next element and we can access next's previous. */
		new_el->BASEREFERENCE(prev) = next_el->BASEREFERENCE(prev);
		next_el->BASEREFERENCE(prev) = new_el;
	} 

	/* Set forward pointers. */
	if (new_el->BASEREFERENCE(prev) == 0) {
		/* There is no previous element. Set the head pointer.*/
		head = new_el;
	}
	else {
		/* There is a previous element, set it's next pointer to new_el. */
		new_el->BASEREFERENCE(prev)->BASEREFERENCE(next) = new_el;
	}

	/* Update list length. */
	listLen++;
}

/**
 * \brief Insert an entire list immediatly after an element in this list.
 *
 * Elements are moved, not copied. Afterwards, the other list is Empty. If
 * prev_el is NULL then the elements are Prepended to the front of the list.
 * If prev_el is not in the list then undefined behaviour results. All
 * elements are inserted into the list at once, so this is an O(1) operation.
 */
template <LINKEDLIST_MULTIELEMENT_TEMPDEF> void LinkedList<LINKEDLIST_MULTIELEMENT_TEMPUSE>::
		AddAfter( Element *prev_el, LinkedList<LINKEDLIST_MULTIELEMENT_TEMPUSE> &dl )
{
	/* Do not bother if dl has no elements. */
	if ( dl.listLen == 0 )
		return;

	/* Set the previous pointer of dl.head to prev_el. We do
	 * this regardless of the state of the list. */
	dl.head->BASEREFERENCE(prev) = prev_el; 

	/* Set forward pointers. */
	if (prev_el == 0) {
		/* There was no prev_el, we are inserting at the head. */
		dl.tail->BASEREFERENCE(next) = head;
		head = dl.head;
	} 
	else {
		/* There was a prev_el, we can access previous next. */
		dl.tail->BASEREFERENCE(next) = prev_el->BASEREFERENCE(next);
		prev_el->BASEREFERENCE(next) = dl.head;
	} 

	/* Set reverse pointers. */
	if (dl.tail->BASEREFERENCE(next) == 0) {
		/* There is no next element. Set the tail pointer. */
		tail = dl.tail;
	}
	else {
		/* There is a next element. Set it's prev pointer. */
		dl.tail->BASEREFERENCE(next)->BASEREFERENCE(prev) = dl.tail;
	}

	/* Update the list length. */
	listLen += dl.listLen;

	/* Empty out dl. */
	dl.head = dl.tail = 0;
	dl.listLen = 0;
}

/**
 * \brief Insert an entire list immediately before an element in this list.
 *
 * Elements are moved, not copied. Afterwards, the other list is Empty. If
 * next_el is NULL then the elements are Appended to the end of the list. If
 * next_el is not in the list then undefined behaviour results. All elements
 * are inserted at once, so this is an O(1) operation.
 */
template <LINKEDLIST_MULTIELEMENT_TEMPDEF> void LinkedList<LINKEDLIST_MULTIELEMENT_TEMPUSE>::
		AddBefore( Element *next_el, LinkedList<LINKEDLIST_MULTIELEMENT_TEMPUSE> &dl )
{
	/* Do not bother if dl has no elements. */
	if ( dl.listLen == 0 )
		return;

	/* Set the next pointer of dl.tail to next_el. We do
	 * this regardless of the state of the list. */
	dl.tail->BASEREFERENCE(next) = next_el; 

	/* Set reverse pointers. */
	if (next_el == 0) {
		/* There is no next elememnt. We are inserting at the tail. */
		dl.head->BASEREFERENCE(prev) = tail;
		tail = dl.tail;
	} 
	else {
		/* There is a next element and we can access next's previous. */
		dl.head->BASEREFERENCE(prev) = next_el->BASEREFERENCE(prev);
		next_el->BASEREFERENCE(prev) = dl.tail;
	} 

	/* Set forward pointers. */
	if (dl.head->BASEREFERENCE(prev) == 0) {
		/* There is no previous element. Set the head pointer.*/
		head = dl.head;
	}
	else {
		/* There is a previous element, set it's next pointer to new_el. */
		dl.head->BASEREFERENCE(prev)->BASEREFERENCE(next) = dl.head;
	}

	/* Update list length. */
	listLen += dl.listLen;

	/* Empty out dl. */
	dl.head = dl.tail = 0;
	dl.listLen = 0;
}


/**
 * \brief Detach an element from the list.
 *
 * The element is not deleted. If the element is not in the list, then
 * undefined behaviour results.
 *
 * \returns The element Detached.
 */
template <LINKEDLIST_MULTIELEMENT_TEMPDEF> Element *LinkedList<LINKEDLIST_MULTIELEMENT_TEMPUSE>::
		Detach(Element *el)
{
	/* Set forward pointers to skip over el. */
	if (el->BASEREFERENCE(prev) == 0) 
		head = el->BASEREFERENCE(next); 
	else {
		el->BASEREFERENCE(prev)->BASEREFERENCE(next) =
				el->BASEREFERENCE(next); 
	}

	/* Set reverse pointers to skip over el. */
	if (el->BASEREFERENCE(next) == 0) 
		tail = el->BASEREFERENCE(prev); 
	else {
		el->BASEREFERENCE(next)->BASEREFERENCE(prev) =
				el->BASEREFERENCE(prev); 
	}

	/* Update List length and return element we Detached. */
	listLen--;
	if(listLen<0)
		listLen=0;
	return el;
}

/**
 * \brief Clear the list by deleting all elements.
 *
 * Each item in the list is deleted. The list is reset to its initial state.
 */
template <LINKEDLIST_MULTIELEMENT_TEMPDEF> void LinkedList<LINKEDLIST_MULTIELEMENT_TEMPUSE>::Empty()
{
	Element *nextToGo = 0, *cur = head;
	
	while (cur != 0)
	{
		nextToGo = cur->BASEREFERENCE(next);
		delete cur;
		cur = nextToGo;
	}
	head = tail = 0;
	listLen = 0;
}

/**
 * \brief Clear the list by forgetting all elements.
 *
 * All elements are Abandoned, not deleted. The list is reset to it's initial
 * state.
 */
template <LINKEDLIST_MULTIELEMENT_TEMPDEF> void LinkedList<LINKEDLIST_MULTIELEMENT_TEMPUSE>::Abandon()
{
	head = tail = 0;
	listLen = 0;
}
