#
# tests the ctypes python importation of libfreckle
#

import sys
from random import choice
from ctypes import *
lib=CDLL("./libfreckle.so")
#lib.getInfo(None)

seq1="".join([choice("ACGT") for num in xrange(1000)])
#seq2="".join([choice("ACGT") for num in xrange(1000)])
seq2=seq1+"".join([a for a in reversed(seq1)])
#seq1="".join([a for a in seq1])

dotstore=lib.makeDotComparison(seq1, seq2, 4, 10, 1, 5)

#lib.DumpDotStore(dotstore)
#sys.exit(0)
#lib.DotStoreImageToString.restype = c_char_p

longest=len(seq2)
print "SEQ1:",len(seq1)
print "SEQ2:",len(seq2)



string = lib.DotStoreImageToString(dotstore,len(seq1),len(seq2),longest,2*4)

print "WTF!"

w = len(seq1)
h = len(seq2)

if h>w:
	ht=longest
	wt=int( (float(w)/float(h))*float(ht)   )
else:
	wt=longest
	ht=int( (float(h)/float(w))*float(wt)   )

print "aha",wt,ht,wt*ht
result=string_at(string,wt*ht)
print "yes!"

print result

#lets negative it
data=''
for pixel in result:
	data+=chr(255-ord(pixel))

from PIL import Image

print "PIL!"
image=Image.fromstring("L", (wt,ht), data)
print "PIL2"
image.save("testimage.png")
print "PIL3"
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

