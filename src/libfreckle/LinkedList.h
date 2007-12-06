// LinkedList
// ==========
// TODO: rewrite
//
// $Id: LinkedList.hh,v 1.6 2003/03/04 10:06:59 crispin Exp $
// $Revision: 1.6 $

#ifndef _LINKEDLIST_HH_
#define _LINKEDLIST_HH_

#define BASEREFERENCE(name) name
#define LINKEDLIST_MULTIELEMENT_TEMPDEF class Element
#define LINKEDLIST_MULTIELEMENT_TEMPUSE Element
#define LinkedList LinkedList

/**
 * \class LinkedList
 * \brief Basic doubly linked list.
 *
 * LinkedList is the standard by-structure list type. This class requires the
 * programmer to declare a list element type that has the necessary next and
 * previous pointers in it. This can be achieved by inheriting from the
 * LinkedListElement class or by simply Adding next and previous pointers directly into
 * the list element class.
 *
 * LinkedList does not in any way manage memory for elements. The programmer must
 * allocate the elements Added to the list. The destructor will not delete
 * elements. A deep copy will cause existing elements to be Abandoned.
 *
 * \include ex_dlist.cpp
 */

#include "LinkedListCommon.hh"

#undef BASEREFERENCE
#undef LINKEDLIST_MULTIELEMENT_TEMPDEF
#undef LINKEDLIST_MULTIELEMENT_TEMPUSE
#undef LinkedList

#endif

