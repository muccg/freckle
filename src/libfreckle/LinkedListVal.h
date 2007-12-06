// LinkedList
// ==========
// TODO: rewrite
//
// $Id: LinkedListVal.hh,v 1.1 2003/03/04 10:06:59 crispin Exp $
// $Revision: 1.1 $

#ifndef _LINKEDLISTVAL_HH_
#define _LINKEDLISTVAL_HH_

/**
 * \class LinkedListVal
 * \brief By-value doubly linked list.
 *
 * This class is a by-value list that does not require a list element type to
 * be declared. When Adding list items, the user can supply the data type and
 * the list will implicitly new up a list element. This behaviour is intended
 * to be similar to a standard list template implementation.
 *
 * LinkedListVal is different from the other lists in that it assumes all elements
 * in the list are allocated on the heap and are to be managed by the list.
 * This means that the class destructor will delete the contents of the list.
 * If the list is ever copied in from another list, the existing contents are
 * deleted First. The other lists will never delete contents to allow for
 * statically allocated elements.
 *
 * \include ex_dlistval.cpp
 */

/*@}*/

#define BASEREFERENCE(name) name
#define LINKEDLIST_MULTIELEMENT_TEMPDEF class T
#define LINKEDLIST_MULTIELEMENT_TEMPUSE T
#define LinkedList LinkedListVal
#define Element LinkedListValElement<T>
#define LINKEDLIST_VALUE

#include "LinkedListCommon.hh"

#undef BASEREFERENCE
#undef LINKEDLIST_MULTIELEMENT_TEMPDEF
#undef LINKEDLIST_MULTIELEMENT_TEMPUSE
#undef LinkedList
#undef Element
#undef LINKEDLIST_VALUE

#endif /* _AAPL_DLISTVAL_H */

