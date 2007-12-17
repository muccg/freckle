from math import *
from Bio import SeqIO
from Bio.Seq import Seq
from numpy import array, zeros, int32, vstack, hstack
from PIL import Image, ImageDraw, ImageFont
from pyfreckle import *
import os.path

class DotPlot:
	"""
	\brief Class encapsulating a DotPlot or a sub part thereof
	\details This is the DotPlot class. Its a python class that helps you very efficiently generate dotplots. It can be used to generate
	a small area of a large dotplot or a complete dotplot graph with axis rulers, titles and annotations.
	To use you basically instantiate it with two lists of files. One list of fasta files for the x axis and one list of fasta files
	for the y axis. Then call methods of the class to generate the graphs. You can also save processed information out and read it
	back in later to save processing time. For instance you can generate an entire dot plot for a very large sequence and save it
	to ram or to disk. Then you can use that dot plot info to render different parts of the dot plot later.
	"""
	def __init__(self, xfiles, yfiles):
		"""
		\brief Creates a DotPlot object using two lists of fasta files as the sequences for the x and y axis.
		
		\param xfiles A list of FASTA formatted files to be used for the x axis.
		\param yfiles A list of FASTA formatted files to be used for the y axis.
		"""
		self.filenames=[xfiles,yfiles]
		
		# sequence bounds is a list of lists for each dimension. It contains the length of each sequence. eg. for
		# a single dimension it may be something like [ [51,23], [45] ]
		# this indicates two files, the first one with two sequences, of length 51 and 23. and a second file of a
		# single sequence of length 45
		self.sequencebounds=[[[len(x.seq) for x in SeqIO.parse(open(file),"fasta")] for file in seqfiles] for seqfiles in self.filenames]
		self.sequenceboundids=[[[x.id for x in SeqIO.parse(open(file),"fasta")] for file in seqfiles] for seqfiles in self.filenames]
		
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
		"""
		\brief return the length of the full sequence specified for dimension.
		\param dimension The dimension we want the size of. 0 for x. 1 for y.
		\return an integer representing the full combined sequence length
		"""
		return self.size[dimension]
		
	def AssembleFullSequence(self, dimension=0):
		"""
		\brief return the sequence of a dimension as one single sequence object.
		\details This concatenates all the sequences of all the files together end to end.
		\warning could be VERY memory hungry with large numbers of sequence or files
		\param dimension The dimension we want the full sequence of. 0 for x. 1 for y.
		\return a Seq() object that is the full combined sequence. The Seq object returned has no id.
		"""
		assert(dimension==0 or dimension==1)
		return reduce(lambda x,y: x+y, [reduce(lambda x,y: x+y, [x.seq for x in SeqIO.parse(open(file),"fasta")]) for file in self.filenames[dimension]])
	
	def CalculateDotStore(self, ktup=8, window=16, minmatch=8, mismatch=0):
		"""
		\brief performs the calculation of the dot plot for the entire sequences
		\details Calls the underlying libfreckle C code to calculate the dot plot for the entire sequence
		\warning can be slow with large sequences. 500000 vs 500000 can take a few hours to compute on 3 gigs of 32bit CPU
		\todo make this method able to calculate sub dot plots.
		\todo implement ktuple size, window and mismatch
		\return Nothing
		"""
		assert(minmatch>=ktup)
		assert(ktup>=4)
		assert(window>=ktup)
		self.fullseq=[self.AssembleFullSequence(i).data for i in [0,1]]
		#self.fullseq[1]=self.fullseq[1][::-1]
		# forward and reverse dot store
		self.dotstore=[makeDotComparison(self.fullseq[0],self.fullseq[1],ktuplesize=ktup,window=2*ktup,minmatch=ktup),makeDotComparison(self.fullseq[0],self.fullseq[1][::-1],ktuplesize=ktup,window=window,minmatch=minmatch,mismatch=mismatch)]
	
	def IndexDotStore(self):
		"""
		\brief Indexes the calculated DotStore
		"""
		# create indexes
		[dots.CreateIndex() for dots in self.dotstore]
	
	def MakeAverageGrid(self,scale):
		"""
		\brief Calculates the reduced score grid
		\details Uses a fairly optimised algorithm to convert the calculated dot store into a grid of averaged values. These
		values represent the count of how many dots would be in that area.
		\todo allow sub windows to be averaged
		\param scale the scale value for the sizing of the grid. >1 to shrink. eg 10 means the final averaged grid
		will be 1/10th the size of the full dotplot.
		\return Nothing
		"""
		assert(self.dotstore)
		self.scale=scale
		grid=[DotGrid(),DotGrid()]
		if len(self.fullseq[0])>len(self.fullseq[1]):
			# wide image
			lenseq1=ceil(scale)*len(self.fullseq[0])/scale
			lenseq2=ceil(scale)*len(self.fullseq[1])/scale
		else:
			# tall image
			lenseq1=ceil(scale)*len(self.fullseq[0])/scale
			lenseq2=ceil(scale)*len(self.fullseq[1])/scale
			
		lenseq1=len(self.fullseq[0])
		lenseq2=len(self.fullseq[1])
			
		[g.Calculate(dots,0,0,lenseq1,lenseq2,scale,10) for g,dots in zip(grid,self.dotstore)]
		grid[1].FlipInplace()
		grid[0].AddInplace(grid[1])
		self.grid=grid[0]
		
	def MakeImage(self):
		"""
		\brief Makes a DotPlot image from the averaged grid data
		"""
		string=self.grid.ToString()
		image=Image.fromstring("L", (self.grid.GetWidth(),self.grid.GetHeight()), string).convert("RGB")
		
		self.DrawBounds(image)
		image=self.AddAxis(image)
		
		return image
	
	def AddAxis(self,image):
		"""
		\brief add the image if an axis and sequence annotations to the dotplot image
		\detail takes the graph image and creates a new image with axis and annotations. the new image will be larger than the old
		\param image the pre calculated image of the dotplot
		\return the annotated image
		"""
		size=image.size
		font = ImageFont.truetype("/usr/share/fonts/truetype/freefont/FreeSans.ttf", 10)
		
		titlefont = ImageFont.truetype("/usr/share/fonts/truetype/freefont/FreeSansBold.ttf", 12)
		namelist=["",""]
		if len(self.filenames[0])==1:
			namelist[0]=os.path.basename(self.filenames[0][0])
		else:
			namelist[0]=", ".join([os.path.basename(a) for a in self.filenames[0][:-1]])
			namelist[0]=namelist[0]+" and "+os.path.basename(self.filenames[0][-1])
		if len(self.filenames[1])==1:
			namelist[1]=os.path.basename(self.filenames[1][0])
		else:
			namelist[1]=", ".join([os.path.basename(a) for a in self.filenames[1][:-1]])
			namelist[1]=namelist[1]+" and "+os.path.basename(self.filenames[1][-1])
		
		title="Dot Plot of " + namelist[0] + " versus "+ namelist[1]
		titlewidth,titleheight=titlefont.getsize(title)
		
		TITLEGAP=10
		
		XOFF,YOFF=10,10						# border
		
		SIDEGUTTER=30
		
		MAJORSCALE=10000
		MINORSCALE=1000
		
		files=[[os.path.basename(a) for a in self.filenames[axis]] for axis in (0,1)]
		fileymaxlen=max([font.getsize(a)[0] for a in files[1]])
		
		
		# find the side width
		ids=[]
		[[ids.append(a) for a in b] for b in self.sequenceboundids[1]]		# 1 is dimension (y)
		idymaxlen=max([font.getsize(a)[0] for a in ids])
		
		# work out the largest sequence number
		lenseq1=len(self.fullseq[0])
		lenseq2=len(self.fullseq[1])
		digits=len(str(max(lenseq1,lenseq2)))			#this is the maximum number of digits in the sequence string
		
		axistextwidth,axistextheight=font.getsize("0"*digits)		#how wide the digit string is in pixels, and how high
		ticklength=10
		
		aisize=(size[0]+axistextwidth+ticklength+2+2*XOFF+idymaxlen+fileymaxlen+2*SIDEGUTTER,size[1]+2+2*YOFF+axistextheight+ticklength+titleheight+TITLEGAP)
		axisimage=Image.new("RGB",aisize,(255,255,255))
		
		dc = ImageDraw.Draw(axisimage,"RGBA")
		dc.rectangle( [(axistextwidth+ticklength+XOFF,axistextheight+ticklength+YOFF+titleheight+TITLEGAP),(aisize[0]-XOFF-1-idymaxlen-fileymaxlen-2*SIDEGUTTER,aisize[1]-YOFF-1)], fill=(0,0,0,255) )
		dc.rectangle( [(axistextwidth+ticklength+XOFF+2,axistextheight+ticklength+YOFF+titleheight+TITLEGAP+2),(aisize[0]-XOFF+1-idymaxlen-fileymaxlen-2*SIDEGUTTER,aisize[1]-YOFF+1)], fill=(0,0,0,255) )
		axisimage.paste(image,(axistextwidth+ticklength+1+XOFF,axistextheight+ticklength+1+YOFF+titleheight+TITLEGAP))
		
		# VERTICAL SCALE
		ticksevery=MAJORSCALE
		scale=self.scale
		value=0.0
		count=0
		while(value<size[1]):
			dc.line( [ (axistextwidth+XOFF,axistextheight+ticklength+value+YOFF+titleheight+TITLEGAP),(axistextwidth+ticklength+XOFF, axistextheight+ticklength+value+YOFF+titleheight+TITLEGAP) ], fill=(0,0,0,255))
			
			text=str(count)
			w,h=font.getsize(text)
			dc.text( (XOFF+axistextwidth-w,value-axistextheight/2+1+YOFF+axistextheight+ticklength+titleheight+TITLEGAP), text, fill=(0,0,0,255), font=font )
			
			value+=float(ticksevery)/scale
			count+=ticksevery
			
		minor=MINORSCALE
		value=0.0
		while(value<size[1]):
			dc.line( [ (axistextwidth+ticklength/2+XOFF,axistextheight+ticklength+value+YOFF+titleheight+TITLEGAP),(axistextwidth+ticklength+XOFF, axistextheight+ticklength+value+YOFF+titleheight+TITLEGAP) ], fill=(0,0,0,255))
			value+=float(minor)/scale
	
		# HORIZONTAL SCALE
		ticksevery=MAJORSCALE
		scale=self.scale
		value=0.0
		count=0
		while(value<size[0]):
			dc.line( [ (value+XOFF+ticklength+axistextwidth,axistextheight+YOFF+titleheight+TITLEGAP),(value+XOFF+ticklength+axistextwidth,axistextheight+ticklength+YOFF+titleheight+TITLEGAP)], fill=(0,0,0,255))
			#dc.line( [ (axistextwidth+XOFF,value+YOFF),(axistextwidth+ticklength+XOFF, value+YOFF) ], fill=(0,0,0,255))
			#dc.text( (XOFF,value-axistextheight/2+1+YOFF), (" "*digits+str(count))[-digits:], fill=(0,0,0,255) )
			text=str(count)
			w,h=font.getsize(text)
			dc.text( (value+XOFF+ticklength+axistextwidth-w/2,YOFF+titleheight+TITLEGAP), text, fill=(0,0,0,255), font=font)
			
			value+=float(ticksevery)/scale
			count+=ticksevery
			
		minor=MINORSCALE
		value=0.0
		while(value<size[0]):
			dc.line( [ (value+XOFF+ticklength+axistextwidth,axistextheight+ticklength/2+YOFF+titleheight+TITLEGAP),(value+XOFF+ticklength+axistextwidth,axistextheight+ticklength+YOFF+titleheight+TITLEGAP)], fill=(0,0,0,255))
			#dc.line( [ (axistextwidth+ticklength/2+XOFF,value+YOFF),(axistextwidth+ticklength+XOFF, value+YOFF) ], fill=(0,0,0,255))
			value+=float(minor)/scale
	
		#title
		dc.text( ((aisize[0]-titlewidth)/2,YOFF), title, fill=(0,0,0,255), font=titlefont)
			
		# size sequence annotations
		axis=1
		files=self.filenames
		boundids=self.sequenceboundids
		seqys=[0]+self.seqybounds
		seqyoff=1
		for file in boundids[axis]:
			for ids in file:
				idy=seqys[seqyoff]-(seqys[seqyoff]-seqys[seqyoff-1])/2
				
				textw,texth=font.getsize(ids)
				x=aisize[0]-XOFF-idymaxlen-fileymaxlen-SIDEGUTTER-SIDEGUTTER/2
				y=YOFF+titleheight+TITLEGAP+ticklength+axistextheight+idy-texth/2
				dc.text( (x,y), ids, fill=(0,0,255,255), font=font)
				
				seqyoff+=1
		
		# file annotations
		axis=1
		files=[os.path.basename(a) for a in self.filenames[axis]]
		fbounds=[0]+[int(a/self.scale) for a in self.globalfilebounds[axis]]
		fyoff=1
		for file in files:
			idy=fbounds[fyoff]-(fbounds[fyoff]-fbounds[fyoff-1])/2
			
			fnamew,fnameh=font.getsize(file)
			x=aisize[0]-XOFF-fileymaxlen
			tops=YOFF+titleheight+TITLEGAP+ticklength+axistextheight
			y=tops+idy-fnameh/2
			dc.text( (x,y), file, fill=(255,0,0,255), font=font)
			
			# draw the lines
			dc.line( [ (aisize[0]-XOFF-fileymaxlen-SIDEGUTTER, fbounds[fyoff-1]+tops+1), (aisize[0]-XOFF-fileymaxlen-2*SIDEGUTTER/3, fbounds[fyoff-1]+tops+1) ], fill=(255,0,0,128) )
			dc.line( [ (aisize[0]-XOFF-fileymaxlen-SIDEGUTTER, fbounds[fyoff]+tops-1), (aisize[0]-XOFF-fileymaxlen-2*SIDEGUTTER/3, fbounds[fyoff]+tops-1) ], fill=(255,0,0,128) )
			dc.line( [ (aisize[0]-XOFF-fileymaxlen-2*SIDEGUTTER/3, fbounds[fyoff-1]+tops+1), (aisize[0]-XOFF-fileymaxlen-2*SIDEGUTTER/3, fbounds[fyoff]+tops-1)], fill=(255,0,0,128))
			dc.line( [ (aisize[0]-XOFF-fileymaxlen-2*SIDEGUTTER/3,y+fnameh/2),(aisize[0]-XOFF-fileymaxlen-SIDEGUTTER/3,y+fnameh/2)],fill=(255,0,0,128))
			
			fyoff+=1
		
		# horizontal annotations at the image bottom
		
		
			
		return axisimage
		
		
	
	def DrawBounds(self,image,filebound=(255,0,0,128),seqbound=(0,0,255,128)):
		"""
		\deprecated
		"""
		dc = ImageDraw.Draw(image,"RGBA")
		
		#dc.text((10, 25), "world", font=font,fill=(0,0,0,255))
		
		xstart=0
		ystart=0
		xend=len(self.fullseq[0])
		yend=len(self.fullseq[1])
		scale=self.scale
		print "SCALE",self.scale
		
		#find sequence bounds and store their relative position to this image
		seqxbounds=[int(float(value-xstart)/scale) for value in reduce(lambda a,b:a+b,self.globalsequencebounds[0]) if value>=xstart and value<=xend]
		seqybounds=[int(float(value-ystart)/scale) for value in reduce(lambda a,b:a+b,self.globalsequencebounds[1]) if value>=ystart and value<=yend]
			
		self.seqybounds=seqybounds[:]
			
		# find bounds that are in this range and store their relative position in this image
		filexbounds=[int(float(value-xstart)/scale) for value in self.globalfilebounds[0] if value>=xstart and value<=xend]
		fileybounds=[int(float(value-ystart)/scale) for value in self.globalfilebounds[1] if value>=ystart and value<=yend]
		
		[seqxbounds.remove(x) for x in filexbounds]
		[seqybounds.remove(y) for y in fileybounds]
					
		for xbound in seqxbounds:
			dc.line( [ (xbound,0), (xbound,image.size[1]) ], fill=seqbound )
		for ybound in seqybounds:
			dc.line( [ (0,ybound), (image.size[0],ybound) ], fill=seqbound )
		
		for xbound in filexbounds:
			dc.line( [ (xbound,0), (xbound,image.size[1]) ], fill=filebound )
		for ybound in fileybounds:
			dc.line( [ (0,ybound), (image.size[0],ybound) ], fill=filebound )
			
		return image
		
	def GetSubSequence(self, dimension, start, end):
		"""
		\brief return an assembled subsequence of the dotplot
		\detail return a sequence that is a combination of one or more sequences, from global seq offset 'start' to 'end'
		\param dimension the dimension to retrieve from. 0 for x. 1 for y
		\param start the start sequence number
		\param end the end sequence number
		\return a new sequence object
		"""
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
