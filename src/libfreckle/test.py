#
# tests the ctypes python importation of libfreckle
#

from ctypes import *
lib=CDLL("./libfreckle.so")
lib.getInfo(None)

#seq1="GATAATGACTGACTGACTGACTGACTGACGACATGGCGCGAGTCGTGATGATATTACTATATATTATATGTGACTGACTGCTGATGCGTATGCGTATGCAGTCTGA"
#seq2="GCTAGCTGTATATATCGCTAGATCTCTTATTATATATAGTGTCGCGTCGTTATGACTATATTATATAAAATACTGCAGTACTGACTGGCGCGCATTCTCTTACTGATGATGCTGCAGTCTGATGCTGACTG"

#dotstore=lib.makeDotComparison(seq1, seq2, 4, 10, 1, 5)

#num=lib.GetNumDots(dotstore)
#for i in xrange(num):
	#print i,lib.GetDotX(dotstore,i),lib.GetDotY(dotstore,i),lib.GetDotLength(dotstore,i)
	
	
#codetable={
	#'G':{
		#'G':{ 'U':'G', 'C':'G', 'A':'G', 'G':'G' },
		#'A':{ 'U':'D', 'C':'D', 'A':'E', 'G':'E' },
		#'C':{ 'U':'A', 'C':'A', 'A':'A', 'G':'A' },
		#'U':{ 'U':'V', 'C':'V', 'A':'V', 'G':'V' }
	#}, 'A':{
		#'G':{ 'U':'S', 'C':'S', 'A':'R', 'G':'R' },
		#'A':{ 'U':'N', 'C':'N', 'A':'K', 'G':'K' },
		#'C':{ 'U':'T', 'C':'T', 'A':'T', 'G':'T' },
		#'U':{ 'U':'I', 'C':'I', 'A':'I', 'G':'M' }
	#}, 'C':{
		#'G':{ 'U':'R', 'C':'R', 'A':'R', 'G':'R' },
		#'A':{ 'U':'H', 'C':'H', 'A':'Q', 'G':'Q' },
		#'C':{ 'U':'P', 'C':'P', 'A':'P', 'G':'P' },
		#'U':{ 'U':'L', 'C':'L', 'A':'L', 'G':'L' }
	#}, 'U':{
		#'G':{ 'U':'C', 'C':'C', 'A':'-', 'G':'W' },
		#'A':{ 'U':'Y', 'C':'Y', 'A':'-', 'G':'-' },
		#'C':{ 'U':'S', 'C':'S', 'A':'S', 'G':'S' },
		#'U':{ 'U':'F', 'C':'F', 'A':'L', 'G':'L' }
	#}
#}

#outs=""
#for p1 in 'ACGU':
	#for p2 in 'ACGU':
		#for p3 in 'ACGU':
			#outs=outs+codetable[p1][p2][p3]
			
#print outs

