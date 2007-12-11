from Bio import SeqIO
from Bio.Seq import Seq
from numpy import array, zeros, int32, vstack, hstack
from PIL import Image, ImageDraw

class DotPlot:
	def __init__(self, xfiles, yfiles):
		self.filenames=[xfiles,yfiles]
		
		# sequence bounds is a list of lists for each dimension. It contains the length of each sequence. eg. for
		# a single dimension it may be something like [ [51,23], [45] ]
		# this indicates two files, the first one with two sequences, of length 51 and 23. and a second file of a
		# single sequence of length 45
		self.sequencebounds=[[[len(x.seq) for x in SeqIO.parse(open(file),"fasta")] for file in seqfiles] for seqfiles in self.filenames]
		
		# the total sequence count of each file. in the above example it would be [ 74, 45 ] (74 = 51+23)
		self.filebounds=[[sum(c) for c in bounds] for bounds in self.sequencebounds ]
				
		# the total size of all the sequences. in the above example the value would be 119. There is one number for each dimension
		self.size=[sum([a for a in fileoff]) for fileoff in self.filebounds]
		
		def accumulate(sequence):
			"""Takes a sequence of values like [10,4,25] and adds each to the next to create the ascending list [10,14,39]"""
			out=[]
			for a in xrange(len(sequence)):
				out.append(sum(sequence[0:a+1]))
			return out
		
		
		#OK all those were relative. Now we are going to create some absolute indexes in 'global' coordinates
		self.globalfilebounds=[ accumulate(a) for a in self.filebounds ]
		self.globalsequencebounds=[ [accumulate(b) for b in a] for a in self.sequencebounds]
		
		# add the file offsets so they truly are global instead of local to the individual file
		self.globalsequencebounds=[ [[b+off for b in a] for a,off in zip(gsb,[0]+gfb[:-1])] for gsb,gfb in zip(self.globalsequencebounds, self.globalfilebounds)]
		
	def GetSequenceLength(self,dimension):
		return self.size[dimension]
		
	def AssembleFullSequence(self, dimension=0):
		"""return the sequence of a dimension as one single sequence object. This concatenates all the sequences of all the files together end to end.
		WARNING: could be VERY memory hungry with large numbers of sequence or files"""
		assert(dimension==0 or dimension==1)
		return reduce(lambda x,y: x+y, [reduce(lambda x,y: x+y, [x.seq for x in SeqIO.parse(open(file),"fasta")]) for file in self.filenames[dimension]])
	
	def DrawDotPlot(self,xstart,ystart,xend=None,yend=None,window=1,filebound=(255,0,0),seqbound=(0,255,0)):
		"""Draw a section of the dot plot including file bounds and sequence bounds as coloured lines"""
		assert(window)
		assert(xstart>=0 and xstart<=self.size[0])
		assert(ystart>=0 and ystart<=self.size[1])
		
		#Python 2.5
		#xend= self.size[0] if xend==None else xend
		#yend= self.size[1] if yend==None else yend
		
		#Python2.4
		xend=(xend==None and self.size[0] or xend)
		yend=(yend==None and self.size[1] or yend)
				
		assert(xend>xstart and xend<=self.size[0])
		assert(yend>ystart and yend<=self.size[1])
		
		dotplot=self.MakeDotPlotSubImage( self.GetSubSequence(0,xstart,xend), self.GetSubSequence(0,ystart,yend), window )
		
		#draw on our file boundary lines
	
		# make an RGB image of the same size
		im=dotplot.convert("RGB")
		
		dc = ImageDraw.Draw(im)
		
		#find sequence bounds and store their relative position to this image
		applicablexbounds=[value-xstart for value in reduce(lambda a,b:a+b,self.globalsequencebounds[0]) if value>=xstart and value<=xend]
		applicableybounds=[value-ystart for value in reduce(lambda a,b:a+b,self.globalsequencebounds[1]) if value>=ystart and value<=yend]
			
		print applicablexbounds
			
		for xbound in applicablexbounds:
			dc.line( [ (xbound,0), (xbound,dotplot.size[1]) ], fill=seqbound )
		for ybound in applicableybounds:
			dc.line( [ (0,ybound), (dotplot.size[0],ybound) ], fill=seqbound )
		
		# find bounds that are in this range and store their relative position in this image
		applicablexbounds=[value-xstart for value in self.globalfilebounds[0] if value>=xstart and value<=xend]
		applicableybounds=[value-ystart for value in self.globalfilebounds[1] if value>=ystart and value<=yend]
		
		
		for xbound in applicablexbounds:
			dc.line( [ (xbound,0), (xbound,dotplot.size[1]) ], fill=filebound )
		for ybound in applicableybounds:
			dc.line( [ (0,ybound), (dotplot.size[0],ybound) ], fill=filebound )
		
		return im
		
		
	def MakeDotPlotSubImage(self, xseq, yseq, windowsize=1):
		"""Here we create a dot plot of the xseq vs the yseq
		Windowsize is how many in a row need to indicate a point
		negative window size looks for reverse sequences
		"""
		assert(windowsize)
		#truthmatrix=255-self.CompareSequences(xseq, yseq)*255
		# first get out matching boolean array
		truthmatrix=self.CompareSequences(xseq, yseq)
		(y,x)=truthmatrix.shape
		
		if windowsize>1:
			#forward match looking
			clip=windowsize-1
			
			# move the window across the grid from top left to bottom right accumulating the result through a bitwise AND
			accum=reduce( lambda x,y: x&y, (truthmatrix[start:yend, start:xend] for start,xend,yend in zip(range(windowsize),range(x-clip,x+2),range(y-clip,y+2))))
		else:
			#reverse match looking
			clip=-windowsize-1
			
			# move the window across the grid from bottom left to top right accumulating the result through a bitwise AND
			accum=reduce( lambda x,y: x&y, (truthmatrix[ystart:yend, xstart:xend] for xstart,ystart,xend,yend in zip(range(-windowsize),reversed(range(-windowsize)),range(x-clip,x+2),reversed(range(y-clip,y+1)))))
		
		# turn it into an image
		imagematrix=255-accum*255
		(y,x)=imagematrix.shape
		
		# convert to a Luminance image
		data=imagematrix.tostring()
		return Image.fromstring("L",(x,y),data)
		
	def GetSubSequence(self, dimension, start, end):
		"""return a sequence that is a combination of one or more sequences, from global seq offset 'start' to 'end'"""
		assert(dimension==0 or dimension==1)		# atm only 2D
		assert(start>=0)
		assert(start<=self.size[dimension])
		assert(end>=0)
		assert(end<=self.size[dimension])
		assert(start<end)
		
		def getoffset(valuelist, value):
			"""returns the SLICE offset of the SLICE POINT of the last value in the ASCENDING list valuelist that is LESS THAN the passed in value. keh?
			if the the passed in values were [10,15,30,47],17 the return value will be 2, because slice position 2 marks the point
			between 15 and 30, where 17 would sit
			"""
			#quick assetion that this is an ascending list
			for a,b in zip(valuelist[:-1],valuelist[1:]):
				assert(a<b)
			
			#TODO: what if value is already in valuelist
				
			return sorted(valuelist[:]+[value]).index(value)
			
		# first we have to work out all the files that are involved. in small files this may be a lot. in large files this will probably be a few
		filenames=self.filenames[dimension]
		
		# the start end end files (as slice offset values)
		startoff,endoff=getoffset(self.globalfilebounds[dimension],start),getoffset(self.globalfilebounds[dimension],end)+1
		outseqlength=end-start
				
		outseq=Seq("")
		for filenum in xrange(startoff, endoff):
			#process file filenum
			seqnum=0
			for seqrecord in SeqIO.parse(open(filenames[filenum]),"fasta"):
				seq=seqrecord.seq
				highval=self.globalsequencebounds[dimension][filenum][seqnum]
				seqsize=self.sequencebounds[dimension][filenum][seqnum]
				if highval>start:
					#this subsequence has begun
					# the start value within this sequence
					startval=seqsize-highval+start
					
					#if we haven't ended yet
					if len(outseq)<outseqlength:
						# if we have a negative start value then the sub sequence started in a previous sequence
						if startval<0:
							#is the end in this sequence?
							relendpos=outseqlength+startval
							if relendpos<seqsize:
								#yes ends in this sequence
								outseq+=seq[:relendpos]
							else:
								#includes the whole sequence
								outseq+=seq
						else:
							#do we also end here?
							if outseqlength+startval<highval:
								#yes. start and end in one
								outseq+=seq[startval:startval+outseqlength]
							else:
								#no include the whole sequence
								outseq+=seq[startval:]
					
					
				seqnum+=1
					
		assert(len(outseq) == end-start)
		
		return outseq






#	
# test suite
#
import unittest, random

class TestDotPlot(unittest.TestCase):
	NUMAXIS=2
	
	def setUp(self):
		from glob import glob
		self.filelist=glob("TestFastaFiles/*.fasta")
	
	def notestBounds(self):
		"""Test the generation and metadata"""
		dp=DotPlot(self.filelist, self.filelist)
		
		# loop over each axis
		for axis in xrange(self.NUMAXIS):
			fb=dp.filebounds[axis]
			gfb=dp.globalfilebounds[axis]
			sb=dp.sequencebounds[axis]
			gsb=dp.globalsequencebounds[axis]
			
			# make sure each gfb entry is the sum of all the previous fb entries
			for num in xrange(len(gfb)):
				self.assertEqual(gfb[num], sum(fb[:num+1]))
				
			# Make sure each sb entry sums to the relevant fb entry
			for num in xrange(len(sb)):
				self.assertEqual(sum(sb[num]),fb[num])
				
			# make sure the final entry in each gsb is equal to the relevant entry from gfb
			for num in xrange(len(gsb)):
				self.assertEqual(gsb[num][-1],gfb[num])
				
			# make sure the accumulation of each sb is the entry in the gsb
			# first flatten the list of lists into a single list
			flatsb,flatgsb=[],[]
			[[flatsb.append(c) for c in b] for b in sb]
			[[flatgsb.append(c) for c in b] for b in gsb]
			for num in xrange(len(flatgsb)):
				self.assertEquals(flatgsb[num], sum(flatsb[:num+1]))
				
	def notestGetSubSequence(self):
		"""Test the retreival of subsequences"""
		dp=DotPlot(self.filelist, self.filelist)
		
		# pick a random subsequence within the range of the files and get that subseq
		# make sure that subsequence is actually the subsequence it should be
		
		seqlen=dp.GetSequenceLength(0)			# total number of sequences
		fullseq=dp.AssembleFullSequence(0)
	
		for testnum in xrange(10):
			startandend=[random.randint(0,seqlen-1) for num in range(2) ]
			start,end=min(startandend),max(startandend)
			
			#WARNING this is extremely memory hungry but is the most robust way of testing
			#we assemble the complete sequence string so we can test subsequence against a slice of the full sequence
			subseqcomparison=fullseq[start:end]
			subseq=dp.GetSubSequence(0,start,end)
			self.assertEquals(subseq.tostring(),subseqcomparison.tostring())
		
	def notestCompareSequences(self):
		"""Test the matrix comparison of sequences"""
		dp=DotPlot(self.filelist, self.filelist)
		
		fullseq=dp.AssembleFullSequence(0)
		seqlen=len(fullseq)
		
		for testnum in xrange(2):
			startandend=[random.randint(0,seqlen-1) for num in range(2) ]
			sx,ex=min(startandend),max(startandend)
			
			startandend=[random.randint(0,seqlen-1) for num in range(2) ]
			sy,ey=min(startandend),max(startandend)
			
			matrix = dp.CompareSequences(fullseq[sx:ex],fullseq[sy:ey])
			
			#make sure its the right shape. shape is returned as height,width (because its [outside [inside] ])
			self.assertEquals( matrix.shape, (ey-sy, ex-sx) )
			
			#make sure all the trues and falses are correct
			for xpos in xrange(ex-sx):
				for ypos in xrange(ey-sy):
					self.assertEquals( matrix[ypos][xpos], fullseq[sx+xpos]==fullseq[sy+ypos] )
					
			
	def notestMakeDotPlotSubImage(self):
		dp=DotPlot(self.filelist, self.filelist)
		
		fullseq=dp.AssembleFullSequence(0)
		seqlen=len(fullseq)
		
		image=dp.MakeDotPlotSubImage(fullseq,fullseq)
		image.save("test.png")
		
	def notestBigPlot(self):
		dp=DotPlot(self.filelist, self.filelist)
		
		image=dp.MakeDotPlotSubImage(dp.GetSubSequence(0,0,1024),dp.GetSubSequence(1,0,512),-5)
		image.save("test.png")
		
	def notestPlotWithBorders(self):
		dp=DotPlot(self.filelist, self.filelist)
		
		image=dp.DrawDotPlot(0,0,window=4)
		image.save("test.png")
		
if __name__ == '__main__':
    unittest.main()
