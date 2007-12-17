from pyfreckle import *
from random import choice, randint
from time import time
from PIL import Image

SIZE=5000.
OUTPUT=1000.
OUTPUTV=1000.

# run a dot comparison
seq1="".join([choice('ACGT') for x in xrange(int(SIZE))])
#seq2="".join([choice('ACGT') for x in xrange(1000)])
seq2=seq1		#+seq1[::-1]
print "seq1 = %d base pairs"%len(seq1)
print "seq2 = %d base pairs"%len(seq2)

start=time()

print "calculating the forward dot plot..."
t=time()
dots=makeDotComparison(seq1,seq2)
print "done in", time()-t, "seconds"
print "%d dot points"%len(dots)

print "calculating the backward dotplot..."
t=time()
dotsbackwards=makeDotComparison(seq1,seq2[::-1])
print "done in", time()-t, "seconds"
print "%d dot points"%len(dotsbackwards)

print "indexing forward store..."
t=time()
dots.CreateIndex()
print "done in", time()-t, "seconds"

print "indexing backward store..."
t=time()
dotsbackwards.CreateIndex()
print "done in", time()-t, "seconds"

# make an averaged grid
print "making averaged grid %dx%d for forward store..."%(OUTPUT,OUTPUTV)
t=time()
grid=DotGrid()
grid.Calculate(dots,0,0,len(seq1),len(seq2),SIZE/OUTPUT,10)
print "done in", time()-t, "seconds"

print "making averaged grid %dx%d for backward store..."%(OUTPUT,OUTPUTV)
t=time()
gridbackwards=DotGrid()
gridbackwards.Calculate(dotsbackwards,0,0,len(seq1),len(seq2),SIZE/OUTPUT,10)
print "done in", time()-t, "seconds"

print "joining grids..."
t=time()
gridbackwards.FlipInplace()
grid.AddInplace(gridbackwards)
print "done in", time()-t, "seconds"
	 
print "calculating histogram..."
t=time()
grid.CalculateHistogram()
print "done in", time()-t, "seconds"
	   
# saving as image
print "saving image.png"
t=time()
string=grid.ToString()
image=Image.fromstring("L", (int(OUTPUT),int(OUTPUTV)), string)
image.save("image.png")
print "done in", time()-t, "seconds"

print "TOTAL TIME:",time()-start,"seconds"

