#
# tests the ctypes python importation of libfreckle
#

from ctypes import *
lib=CDLL("./libfreckle.so")
#lib.getInfo(None)

seq1="GATAATGACTGACTGACTGACTGACTGACGACATGGCGCGAGTCGTGATGATATTACTATATATTATATGTGACTGACTGCTGATGCGTATGCGTATGCAGTCTGA"
seq2="GCTAGCTGTATATATCGCTAGATCTCTTATTATATATAGTGTCGCGTCGTTATGACTATATTATATAAAATACTGCAGTACTGACTGGCGCGCATTCTCTTACTGATGATGCTGCAGTCTGATGCTGACTG"

dotstore=lib.makeDotComparison(seq1, seq2, 4, 10, 1, 5)

num=lib.GetNumDots(dotstore)
for i in xrange(num):
	print i,lib.GetDotX(dotstore,i),lib.GetDotY(dotstore,i),lib.GetDotLength(dotstore,i)