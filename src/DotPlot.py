from math import *
from Bio import SeqIO
from Bio.Seq import Seq
from numpy import array, zeros, int32, vstack, hstack
from PIL import Image, ImageDraw, ImageFont, ImageChops
from pyfreckle import *
import os.path
from time import time
from struct import pack, unpack, calcsize
import pickle
from collections import defaultdict

def locatefile(filename):
	"""use the 'locate' shell command to find all the locations of a file"""
	import popen2
	
	stdout, stdin = popen2.popen2("locate %s"%filename, mode='r')
	stdin.close()
	res=[(a[-1]=='\n') and a[:-1] or a for a in stdout.readlines()]
	stdout.close()
	return res
	
def getfont(filename,size):
	"""try really hard to find and load a particular file"""
	try:
		font = ImageFont.truetype(filename, size)
	except IOError, e:
		#find the font file
		basename=os.path.basename(filename)
		locations=locatefile(basename)
		if len(locations)==0:
			raise Exception, "Cannot find %s font"%basename
		
		font=None
		index=0
		while font==None:
			try:
				font = ImageFont.truetype(locations[index],size)
			except IOError, e:
				index+=1
				if index>=len(locations):
					raise Exception, "Cannot find %s font"%basename
	return font
	
class DotPlotFileError(Exception):
	pass

def decodeseq(seq):
	return ''.join([(a in 'ACGT') and a or 'N'  for a in seq.data.upper()])

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
	def __init__(self, xfiles, yfiles, ktup=8, window=16, minmatch=8, mismatch=0):
		"""
		\brief Creates a DotPlot object using two lists of fasta files as the sequences for the x and y axis.
		
		\param xfiles A list of FASTA formatted files to be used for the x axis.
		\param yfiles A list of FASTA formatted files to be used for the y axis.
		\param ktup The ktuple size for tokenisation
		\param window The size of the window for mismatch comparison
		\param minmatch The minimum length of a match to be included
		\param mismatch The number of mismatching sequences that can be included
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
		
		# our libfreckle calculation tables. One for each dimension
		self.tables=[{},{}]			# one for each dimension. The {} keys are (start,end) pairs with the values being the C tables
		
		# where we put our dotstores
		# indexes is tuple of form (dimension,start,end,compstart,compend)
		# dimension is the classes dimenstion mapped to the dotstores "y" value. This will usually be 1 unless its a special case
		self.dotstore={} 
		
		# the same kind of index but for a calculated average grid
		self.grid={}
		
		# scale is initially unknown
		self.scale=None
		
		# do we autocalc window size
		if window==None:
			window=sum(self.size)/1000
			if window>3:
				window=int(log(window)*3.5)
			else:
				window=8
		
		# save our parameters
		assert(minmatch>=ktup)
		assert(ktup>=4)
		if window<ktup:
			window=ktup
		self.ktup, self.window, self.minmatch, self.mismatch = ktup,window,minmatch,mismatch
		print "Ktuple-size:%d  Window:%d  Minimum-match:%d  Mismatch:%d"%(ktup,window,minmatch,mismatch)
		
		# some expressions to process bounds for a dimension
		self.ProcStart = lambda st: ((st==None) and [0] or [st])[0]
		self.ProcEnd = lambda dim,en: ((en==None) and [self.GetSequenceLength(dim)] or [en])[0]
		
	def Filter(self, length):
		"""
		\brief filter all the stored dotstores to only include matches of at least length
		\param length The minimum allowable match length
		"""
		newds={}
		for key in self.dotstore.keys():
			newds[key]=(self.dotstore[key][0].Filter(length),self.dotstore[key][1].Filter(length))
			
		self.dotstore=newds
		
	def Interpolate(self):
		"""Interpolate all the dot stores"""
		for key in self.dotstore.keys():
			self.dotstore[key][0].Interpolate(self.window)
			self.dotstore[key][1].Interpolate(self.window)
			
		
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
	
	def CreateTables(self, dimension=0, start=None, end=None):
		"""
		\brief creates a set of look up tables for the sequence for the dimension specified
		\details Calls the libfreckle code to make the "C" and "D" tables for the sequence. Stores the result internally
		as state and also returns the object for functional use
		\param dimension the dimension of the sequence to tokenise
		\param start the start index number
		\param end the end index number
		\return the C table pointers in a tuple as (start, end, seqstring, tables)
		"""
		assert(dimension==0 or dimension==1)
		
		start=self.ProcStart(start)
		end=self.ProcEnd(dimension,end)
		
		# munge anything thats not 'ACGT' into a '.'
		subseq=decodeseq(self.GetSubSequence(dimension,start,end))
		
		table=(start,end,subseq,buildMappingTables(subseq, self.ktup))
		self.tables[dimension][(start,end)]=table
		
		return table

	def CreateConservedStore(self, dimension=1):
		"""
		\brief setup the dotplot for conserved region plotting
		\details creates a conserved region dotstore in which to accumulate the conserved areas with extra FASTA files
		"""
		assert len(self.dotstore), "There must be at least one dotstore already calculated for the dotplot prior to creating the conserved store"

		store=self.conservedstore=[DotStore(),DotStore()]

		for s in store:
			s.SetMaxX(self.ProcEnd(dimension,None))
			s.SetMaxY(self.ProcEnd(1-dimension,None))

		self.seqstore=[str(self.AssembleFullSequence(dim)).upper() for dim in (0,1)]

		return store

	def ProcessConservedRegions(self, filename):
		# 1. run the other two dotplots
		#print "calculating conserved dotplots..."
		dp1 = LBDotPlot( self.filenames[0], [filename], ktup=11, window=19, minmatch=11, mismatch=1 )
		dp1.CalculateDotStore()
		
		dp2 = LBDotPlot( self.filenames[1], [filename], ktup=11, window=17, minmatch=11, mismatch=1 )
		dp2.CalculateDotStore()
		#print "done."
		
		conserved = [DotStore(),DotStore()]
		
		#shortcuts for our dotplot results.
		pos,neg = self.dotstore.items()[0][1]
		pos1,neg1=dp1.dotstore.items()[0][1]
		pos2,neg2=dp2.dotstore.items()[0][1]
		
		# lengths
		len0 = len(self.seqstore[0])
		len1 = len(self.seqstore[1])
		len2 = len(dp1.AssembleFullSequence(1))
		
		neg.FlipY(len1)
		#sys.exit()
		
		#print "LEN",len0,len1,len2
		#sys.exit()
		
		#dump=neg
		##print dump, len(dump)
		#print "-----------------------"
		#print "[neg strand]"
		#print "%d dots"%len(dump)
		#for n in range(len(dump)):
			#d = dump[n]
			#print "%d\t%d\t%d"%(d.x,d.y,d.length)
		#sys.exit()
		
		#print "POSITIVE vs POSITIVE/POSITIVE"
		self.Step2( (pos,pos1,pos2), (1,1,1), (len0,len1,len2), conserved[0] )
		
		#print "POSITIVE vs NEGATIVE/NEGATIVE"
		self.Step2( (pos,neg1,neg2), (1,-1,-1), (len0,len1,len2), conserved[0] )
		
		#print "NEGATIVE vs POSITIVE/NEGATIVE"
		self.Step2( (neg,pos1,neg2), (-1,1,-1), (len0,len1,len2), conserved[1] )
		
		#print "NEGATIVE vs NEGATIVE/POSITIVE"
		self.Step2( (neg,neg1,pos2), (-1,-1,1), (len0,len1,len2), conserved[1] )
		
		#print "--------------------------------------"
		
		conserved[1].FlipY(len1)
		neg.FlipY(len1)
		
		return conserved
		
	def CalculateConservedIntersection(self, conservedlist):
		#print [(a,len(a)) for a in self.conservedstore]
		
		conservedstore = self.conservedstore
		X = len(self.seqstore[0])
		Y = len(self.seqstore[1])
		
		# for each plot in the conserved list
		diags = [defaultdict(list), defaultdict(list)]
		
		for cl in conservedlist:
			
			for dpi in (0,1):
				dp = cl[dpi]
				
				# for each diagonal line (defined by x-y)
				for di in range(len(dp)):
					dot = dp[di]
					x,y,l = dot.x, dot.y, dot.length
					
					# set n to be the position up the diagonal
					if x<=y:
						n=x
					else:
						n=y
					
					# add leading edge
					diags[dpi][(x-y)].append( (n,True) )				# add it to the list, or make the list if x-y is not a valid key
				
					# add trailing edge
					diags[dpi][(x-y)].append( (n+l,False) )				# add it to the list, or make the list if x-y is not a valid key
			
			
		# do for both forward and reverse
		for dpi in (0,1):
				
			# at this stage we have diags[V] represent the superposition of the Vth positive diagonal
			# now we can go through each diagonal from -Y to +X and sort over these "edges"
			# we increase a count on a leading edge and decrease a count on a trailing edge, as we move through the edges in order
			# when the count hits len(conservedlist), then its the begining of a point in the intersection
			# when the count drops from len(conservedlist), then its the end of an intersection point
			output={}
			for key in range(-Y,X+1):
				diag = diags[dpi][key]
				if len(diag):
					diag.sort()
					count=0
					start=-1
					end=-1
					parts = []
					for value,edge in diag:
						if edge:
							count+=1
						else:
							count-=1
						
						if count==len(conservedlist):
							# weve hit the number
							start=value
						else:
							# we aren't on the number
							# did we fall from the number?
							if start!=-1:
								end=value
								parts.append( (start,end) )	
								start=-1
							
					output[key]=parts
					
				
			# now we reduce these representations of the diagonal into their cartesian X/Y form
			for key in output.keys():
				# set the starting point of this diagonal
				if key<0:
					x=0
					y=-key
				else:
					x=key
					y=0
					
				# go through starts and ends
				for start,end in output[key]:
					# turn into cartesian dot position
					xpos,ypos,length = x+start, y+start, end-start
					conservedstore[dpi].AddDot(xpos, ypos, length)
					
				
				
			
				
				#conservedstore[0].AddDot(dot.x,dot.y,dot.length)
				
			#for i in range(len(neg)):
				#dot = neg[i]
				#conservedstore[1].AddDot(dot.x,dot.y,dot.length)
			
		
			
		self.dotstore["conserved"]=self.conservedstore
		
	def Step2(self,dotstores, directions, lengths, conservedstore):
		# dotstores should (ds1, ds2, ds3)
		assert len(dotstores)==3
		# directions should be ( +/- 1, +/- 1, +/- 1)
		assert len(directions)==3
		assert False not in [directions[a]==1 or directions[a]==-1 for a in (0,1,2)]
		# lengths should be the lengths of all the sequences
		assert len(lengths)==3
		
		
		ds1,ds2,ds3=dotstores
		dir1,dir2,dir3=directions
		len1,len2,len3=lengths
		
		dottemp={}
		
		for doti in range(len(ds1)):
			#print "%d/%d"%(doti,len(ds1))
			incdots1=[]
			incdots2=[]
			
			dot=ds1[doti]
			x,y,length=dot.x,dot.y,dot.length
			
			#y = (len1-y) if dir1==-1 else y
			#seqx = self.seqstore[0][dot.x:dot.x+dot.length]
			
			x1,x2=x,x+length
			
			# loop over ds2 dots looking for any dots that fall in this region
			for di in range(len(ds2)):
				#print di,"/",len(pos1)
				d=ds2[di]
				if dir2==-1:
					d.y = (len2-d.y)
				dx,dend=d.x,d.x+d.length
				if	not ( (dx<x1 and dend<=x1) or (dx>=x2 and dend>x2) ):
					# we overlap and should be included
					incdots1.append((d.x,d.y,d.length))
			
			#seqy = self.seqstore[1][dot.y:dot.y+dot.length]
			y1,y2=y,y+length
			
			# loop over dp2 dots looking for any positive dots that fall in this region
			for di in range(len(ds3)):
				#print di,"/",len(pos2)
				d=ds3[di]
				if dir3==-1:
					d.y = (len3-d.y)
				dx,dend=d.x,d.x+d.length
				if	not ( (dx<y1 and dend<=y1) or (dx>=y2 and dend>y2) ):
					# we overlap and should be included
					incdots2.append((d.x,d.y,d.length))
			
			#print "INC",incdots1,incdots2
			#if len(incdots1) and len(incdots2):
				#print "FOUND AddDot",incdots1,incdots2
			#print "done."
			
			#print "step 3..."
			
			overlap = lambda a0,a1,b0,b1: not ( ((a0+a1) < b0 ) or ( (b0+b1) < a0 ) ) 
			
			def CalculateOverlap( astart, alength, bstart, blength ):
				aend = astart + alength
				bend = bstart + blength
				
				if bstart >= astart and bstart <= aend:
					if bend < aend:
						return bstart, bend
					return bstart, aend
				
				if bend >= astart and bend <= aend:
					return astart, bend
				
				if bstart < astart and bend > aend:
					return astart, aend
				
				print astart,aend,bstart,bend
				assert False
			
			X,Y,L=0,1,2
			
			if dir2==1 and dir3==1:
				compares = [(a,b) for a in incdots1 for b in incdots2 if overlap(a[Y],a[L],b[Y],b[L])]
			elif dir2==-1 and dir3==1:
				compares = [(a,b) for a in incdots1 for b in incdots2 if overlap(len2-a[Y],a[L],b[Y],b[L])]
			elif dir2==1 and dir3==-1:
				compares = [(a,b) for a in incdots1 for b in incdots2 if overlap(a[Y],a[L],len3-b[Y],b[L])]
			elif dir2==-1 and dir3==-1:
				compares = [(a,b) for a in incdots1 for b in incdots2 if overlap(len2-a[Y],a[L],len3-b[Y],b[L])]
			else:
				assert False, "Unknown direction parameter"
				
			#if len(compares):
				#print "COMPARES",compares
				
			for d1, d2 in compares:
				# compare d1 and d2 on the j axis
				#print d1,d2,(dot.x,dot.y,dot.length)
				
				# cast these to the axis
				xstart, xend = CalculateOverlap( dot.x, dot.length, d1[0], d1[2] )		# [0]=x, [2]=l
				ystart, yend = CalculateOverlap( dot.y, dot.length, d2[0], d2[2] )		# [0]=x, [2]=l
				
				# then cast these to the original "P" dot match
				xoff=xstart-dot.x
				yoff=ystart-dot.y
				
				if xoff > yoff:
					matchstart=xoff
					matchlength=xend-xstart
				elif yoff > xoff:
					matchstart=yoff
					matchlength=yend-ystart
				else:
					#print xoff, yoff,":",xstart,xend,xend-xstart,"|",ystart,yend,yend-ystart
					assert xoff==yoff
					if xend-xstart < yend-ystart:
						matchstart=xoff
						matchlength=xend-xstart
					elif xend-xstart > yend-ystart:
						matchstart=yoff
						matchlength=yend-ystart
					else:
						assert xend-xstart == yend-ystart
						matchstart=yoff
						matchlength=yend-ystart
					
				xs = dot.x+matchstart
				ys = dot.y+matchstart
					
				# add this dot to the conserved dotstore
				#print xend-xstart,yend-ystart,xstart,xend,ystart,yend,"=>",xs,ys,matchlength
				
				x,y,l=xstart,ystart,xend-xstart
				if dottemp.has_key((x,y)):
					#is this longer?
					if l > dottemp[(x,y)]:
						dottemp[(x,y)]=l
				else:
					dottemp[(x,y)]=l
					
				#print "final dot:",x,y,l
					
		# add the dots for real
		for x,y in dottemp.keys():
			conservedstore.AddDot(x,y,dottemp[(x,y)])
		
		

	def DeprecatedProcessConservedRegions(self, filename):
		# 1. run the other two dotplots
		#print "calculating conserved dotplots..."
		dp1 = LBDotPlot( self.filenames[0], [filename] )
		dp1.CalculateDotStore()
		
		dp2 = LBDotPlot( self.filenames[1], [filename] )
		dp2.CalculateDotStore()
		#print "done."
		
		#print "step 2..."
		pos,neg = self.dotstore.items()[0][1]
		pos1,neg1=dp1.dotstore.items()[0][1]
		pos2,neg2=dp2.dotstore.items()[0][1]
		
		##
		## POSITIVE
		##
		print "POS AddDot ----------------------------------------------"
		dottemp={}
			
		for doti in range(len(pos)):
			incdots1=[]
			incdots2=[]
			
			dot=pos[doti]
			x,y,length=dot.x,dot.y,dot.length
			seqx = self.seqstore[0][dot.x:dot.x+dot.length]
			
			x1,x2=x,x+length
			
			# loop over dp1 dots looking for any positive dots that fall in this region
			#print x,y,length
			for di in range(len(pos1)):
				#print di,"/",len(pos1)
				d=pos1[di]
				dx,dend=d.x,d.x+d.length
				if	not ( (dx<x1 and dend<=x1) or (dx>=x2 and dend>x2) ):
					# we overlap and should be included
					incdots1.append((d.x,d.y,d.length))
			
			seqy = self.seqstore[1][dot.y:dot.y+dot.length]
			y1,y2=y,y+length
			
			# loop over dp2 dots looking for any positive dots that fall in this region
			for di in range(len(pos2)):
				#print di,"/",len(pos2)
				d=pos2[di]
				dx,dend=d.x,d.x+d.length
				if	not ( (dx<y1 and dend<=y1) or (dx>=y2 and dend>y2) ):
					# we overlap and should be included
					incdots2.append((d.x,d.y,d.length))
			
			#print "INC",incdots1,incdots2
			if len(incdots1) and len(incdots2):
				print "FOUND AddDot"
			#print "done."
			
			#print "step 3..."
			
			overlap = lambda a,b: not ( ((a[1]+a[2]) < b[1] ) or ( (b[1]+b[2]) < a[1] ) ) 
			
			def CalculateOverlap( astart, alength, bstart, blength ):
				aend = astart + alength
				bend = bstart + blength
				
				if bstart >= astart and bstart <= aend:
					if bend < aend:
						return bstart, bend
					return bstart, aend
				
				if bend >= astart and bend <= aend:
					return astart, bend
				
				if bstart < astart and bend > aend:
					return astart, aend
				
				print astart,aend,bstart,bend
				assert False
					
			
			compares = [(x,y) for x in incdots1 for y in incdots2 if overlap(x,y)]
			for d1, d2 in compares:
				# compare d1 and d2 on the j axis
				#print d1,d2,(dot.x,dot.y,dot.length)
				
				# cast these to the axis
				xstart, xend = CalculateOverlap( dot.x, dot.length, d1[0], d1[2] )
				ystart, yend = CalculateOverlap( dot.y, dot.length, d2[0], d2[2] )
				
				# then cast these to the original "P" dot match
				xoff=xstart-dot.x
				yoff=ystart-dot.y
				
				if xoff > yoff:
					matchstart=xoff
					matchlength=xend-xstart
				elif yoff > xoff:
					matchstart=yoff
					matchlength=yend-ystart
				else:
					#print xoff, yoff,":",xstart,xend,xend-xstart,"|",ystart,yend,yend-ystart
					assert xoff==yoff
					if xend-xstart < yend-ystart:
						matchstart=xoff
						matchlength=xend-xstart
					elif xend-xstart > yend-ystart:
						matchstart=yoff
						matchlength=yend-ystart
					else:
						assert xend-xstart == yend-ystart
						matchstart=yoff
						matchlength=yend-ystart
					
				xs = dot.x+matchstart
				ys = dot.y+matchstart
					
				# add this dot to the conserved dotstore
				#print xend-xstart,yend-ystart,xstart,xend,ystart,yend,"=>",xs,ys,matchlength
				
				x,y,l=xstart,ystart,xend-xstart
				if dottemp.has_key((x,y)):
					#is this longer?
					if l > dottemp[(x,y)]:
						dottemp[(x,y)]=l
				else:
					dottemp[(x,y)]=l
				
		# add the dots for real
		for x,y in dottemp.keys():
			self.conservedstore[0].AddDot(x,y,dottemp[(x,y)])
				
			
				
			print "LEN=",len(compares)

		##
		## NEGATIVE
		##
		print "Neg AddDot ----------------------------------------------"
	
		for doti in range(len(neg)):
			incdots1=[]
			incdots2=[]
			
			dot=neg[doti]
			x,y,length=dot.x,dot.y,dot.length
			seqx = self.seqstore[1][dot.y:dot.y+dot.length]
			
			x1,x2=x,x+length
			
			# loop over dp1 dots looking for any positive dots that fall in this region
			#print x,y,length
			for di in range(len(neg1)):
				#print di,"/",len(pos1)
				d=neg1[di]
				dx,dend=d.x,d.x+d.length
				#print "(%d, %d) : (%d, %d)"%(dx,dend,x1,x2)
				if	not ( (dx<x1 and dend<=x1) or (dx>=x2 and dend>x2) ):
					# we overlap and should be included
					incdots1.append((d.x,d.y,d.length))
			
			print "incdots1 =",incdots1
			
			seqy = self.seqstore[0][dot.x:dot.x+dot.length]
			y1,y2=y,y+length
			
			# loop over dp2 dots looking for any positive dots that fall in this region
			for di in range(len(neg2)):
				#print "-ve:",di,"/",len(neg2)
				d=neg2[di]
				dx,dend=d.x,d.x+d.length
				if	not ( (dx<y1 and dend<=y1) or (dx>=y2 and dend>y2) ):
					# we overlap and should be included
					incdots2.append((d.x,d.y,d.length))
			
			print "incdots2 =",incdots2
			
			if len(incdots1) and len(incdots2):
				print "FOUND AddDot"
			
			#print "done."
			
			#print "step 3..."
			
			overlap = lambda a,b: not ( ((a[1]+a[2]) < b[1] ) or ( (b[1]+b[2]) < a[1] ) ) 
			
			def CalculateOverlap( astart, alength, bstart, blength ):
				aend = astart + alength
				bend = bstart + blength
				
				if bstart >= astart and bstart <= aend:
					if bend < aend:
						return bstart, bend
					return bstart, aend
				
				if bend >= astart and bend <= aend:
					return astart, bend
				
				if bstart < astart and bend > aend:
					return astart, aend
				
				print astart,aend,bstart,bend
				assert False
					
			
			compares = [(x,y) for x in incdots1 for y in incdots2 if overlap(x,y)]
			for d1, d2 in compares:
				# compare d1 and d2 on the j axis
				print "NEG comp:",d1,d2,(dot.x,dot.y,dot.length)
				
				# cast these to the axis
				xstart, xend = CalculateOverlap( dot.x, dot.length, d1[0], d1[2] )
				ystart, yend = CalculateOverlap( dot.y, dot.length, d2[0], d2[2] )
				
				# then cast these to the original "P" dot match
				xoff=xstart-dot.x
				yoff=ystart-dot.y
				
				if xoff > yoff:
					matchstart=xoff
					matchlength=xend-xstart
				elif yoff > xoff:
					matchstart=yoff
					matchlength=yend-ystart
				else:
					assert xoff==yoff
					assert xend-xstart == yend-ystart
					matchstart=yoff
					matchlength=yend-ystart
					
				xs = dot.x+matchstart
				ys = dot.y+matchstart
					
				# add this dot to the conserved dotstore
				#print xend-xstart,yend-ystart,xstart,xend,ystart,yend,"=>",xs,ys,matchlength
				conservedstore.AddDot(xstart,ystart,xend-xstart)
				
			#print "LEN=",len(compares)
			
			
		##
		## Neg vs Pos
		##
		print "Neg vs Pos AddDot ----------------------------------------------"

		for doti in range(len(neg)):
			incdots1=[]
			incdots2=[]
			
			dot=neg[doti]
			x,y,length=dot.x,dot.y,dot.length
			seqx = self.seqstore[1][dot.y:dot.y+dot.length]
			
			x1,x2=x,x+length
			
			# loop over dp1 dots looking for any positive dots that fall in this region
			#print x,y,length
			for di in range(len(neg1)):
				#print di,"/",len(pos1)
				d=neg1[di]
				dx,dend=d.x,d.x+d.length
				#print "(%d, %d) : (%d, %d)"%(dx,dend,x1,x2)
				if	not ( (dx<x1 and dend<=x1) or (dx>=x2 and dend>x2) ):
					# we overlap and should be included
					incdots1.append((d.x,d.y,d.length))
			
			print "incdots1 =",incdots1
			
			seqy = self.seqstore[1][dot.y:dot.y+dot.length]
			y1,y2=y,y+length
			
			# loop over dp2 dots looking for any positive dots that fall in this region
			for di in range(len(pos2)):
				#print di,"/",len(pos2)
				d=pos2[di]
				dx,dend=d.x,d.x+d.length
				if	not ( (dx<y1 and dend<=y1) or (dx>=y2 and dend>y2) ):
					# we overlap and should be included
					incdots2.append((d.x,d.y,d.length))
					
			print "incdots2 =",incdots2
			
			if len(incdots1) and len(incdots2):
				print "FOUND AddDot"
			
			#print "done."
			
			#print "step 3..."
			
			overlap = lambda a,b: not ( ((a[1]+a[2]) < b[1] ) or ( (b[1]+b[2]) < a[1] ) ) 
			
			def CalculateOverlap( astart, alength, bstart, blength ):
				aend = astart + alength
				bend = bstart + blength
				
				if bstart >= astart and bstart <= aend:
					if bend < aend:
						return bstart, bend
					return bstart, aend
				
				if bend >= astart and bend <= aend:
					return astart, bend
				
				if bstart < astart and bend > aend:
					return astart, aend
				
				print astart,aend,bstart,bend
				assert False
					
			
			compares = [(x,y) for x in incdots1 for y in incdots2 if overlap(x,y)]
			for d1, d2 in compares:
				# compare d1 and d2 on the j axis
				print "AddDot NEG/POS comp:",d1,d2,(dot.x,dot.y,dot.length)
				
				# cast these to the axis
				xstart, xend = CalculateOverlap( dot.x, dot.length, d1[0], d1[2] )
				ystart, yend = CalculateOverlap( dot.y, dot.length, d2[0], d2[2] )
				
				# then cast these to the original "P" dot match
				xoff=xstart-dot.x
				yoff=ystart-dot.y
				
				if xoff > yoff:
					matchstart=xoff
					matchlength=xend-xstart
				elif yoff > xoff:
					matchstart=yoff
					matchlength=yend-ystart
				else:
					assert xoff==yoff
					assert xend-xstart == yend-ystart
					matchstart=yoff
					matchlength=yend-ystart
					
				xs = dot.x+matchstart
				ys = dot.y+matchstart
					
				# add this dot to the conserved dotstore
				#print xend-xstart,yend-ystart,xstart,xend,ystart,yend,"=>",xs,ys,matchlength
				self.conservedstore[1].AddDot(xstart,ystart,xend-xstart)
				
			#print "LEN=",len(compares)
			
			
		##
		## Pos vs Neg
		##
		print "Pos vs Neg AddDot ----------------------------------------------"

		for doti in range(len(neg)):
			incdots1=[]
			incdots2=[]
			
			dot=neg[doti]
			x,y,length=dot.x,dot.y,dot.length
			seqx = self.seqstore[1][dot.y:dot.y+dot.length]
			
			x1,x2=x,x+length
			
			# loop over dp1 dots looking for any positive dots that fall in this region
			#print x,y,length
			for di in range(len(pos1)):
				#print di,"/",len(pos1)
				d=pos1[di]
				dx,dend=d.x,d.x+d.length
				if	not ( (dx<x1 and dend<=x1) or (dx>=x2 and dend>x2) ):
					# we overlap and should be included
					incdots1.append((d.x,d.y,d.length))
			
			print "incdots1 =",incdots1
			
			seqy = self.seqstore[0][dot.x:dot.x+dot.length]
			y1,y2=y,y+length
			
			# loop over dp2 dots looking for any positive dots that fall in this region
			for di in range(len(neg2)):
				#print "-ve:",di,"/",len(neg2)
				d=neg2[di]
				dx,dend=d.x,d.x+d.length
				if	not ( (dx<y1 and dend<=y1) or (dx>=y2 and dend>y2) ):
					# we overlap and should be included
					incdots2.append((d.x,d.y,d.length))
			
			print "incdots2 =",incdots2
			
			if len(incdots1) and len(incdots2):
				print "FOUND AddDot"
			
			#print "done."
			
			#print "step 3..."
			
			overlap = lambda a,b: not ( ((a[1]+a[2]) < b[1] ) or ( (b[1]+b[2]) < a[1] ) ) 
			
			def CalculateOverlap( astart, alength, bstart, blength ):
				aend = astart + alength
				bend = bstart + blength
				
				if bstart >= astart and bstart <= aend:
					if bend < aend:
						return bstart, bend
					return bstart, aend
				
				if bend >= astart and bend <= aend:
					return astart, bend
				
				if bstart < astart and bend > aend:
					return astart, aend
				
				print astart,aend,bstart,bend
				assert False
					
			
			compares = [(x,y) for x in incdots1 for y in incdots2 if overlap(x,y)]
			for d1, d2 in compares:
				# compare d1 and d2 on the j axis
				print "NEG comp:",d1,d2,(dot.x,dot.y,dot.length)
				
				# cast these to the axis
				xstart, xend = CalculateOverlap( dot.x, dot.length, d1[0], d1[2] )
				ystart, yend = CalculateOverlap( dot.y, dot.length, d2[0], d2[2] )
				
				# then cast these to the original "P" dot match
				xoff=xstart-dot.x
				yoff=ystart-dot.y
				
				if xoff > yoff:
					matchstart=xoff
					matchlength=xend-xstart
				elif yoff > xoff:
					matchstart=yoff
					matchlength=yend-ystart
				else:
					assert xoff==yoff
					assert xend-xstart == yend-ystart
					matchstart=yoff
					matchlength=yend-ystart
					
				xs = dot.x+matchstart
				ys = dot.y+matchstart
					
				# add this dot to the conserved dotstore
				#print xend-xstart,yend-ystart,xstart,xend,ystart,yend,"=>",xs,ys,matchlength
				self.conservedstore[1].AddDot(xstart,ystart,xend-xstart)
				
			#print "LEN=",len(compares)
			
			
		
		self.dotstore["conserved"]=self.conservedstore

	def DeprecatedProcessConservedRegions(self, filename, dimension=1):
		"""
		\brief process the passed in fasta file and look for conserved areas in relation to this dotplot
		\details loads in the specified fasta file. Then indexes it. Then for each "dot" in this dotplot, find matching regions in the other area. Add the
		matching dot pos and length into the conserved store
		"""
		# get the sequence we are searching for conserved regions
		sequence = str(reduce(lambda x,y: x+y, [x.seq for x in SeqIO.parse(open(filename),"fasta")])).upper()

		# we want to index this 3rd sequence
		tables=buildMappingTables(sequence, self.ktup)			# our tables var = <ctypes.LP_c_void object at 0x??????>

		# our dotstores
		forward, backward = self.dotstore.items()[0][1]

		# for each dot&length in the dotstore
		for dotnum in range(len(forward)):
			print dotnum,"/",len(forward)
			dot = forward.GetDot(dotnum)
			seqx = self.seqstore[0][dot.x:dot.x+dot.length]
			seqy = self.seqstore[1][dot.y:dot.y+dot.length]
			#seqx = self.GetSubSequence(0,dot.x,dot.x+dot.length)
			#seqy = self.GetSubSequence(1,dot.y,dot.y+dot.length)
			assert seqx==seqy

			# now using the tables indexing, search for the longest match from this plot sequence.
			print "findLongestMatch(): sequence=%d seqx=%d x=%d y=%d length=%d"%(len(sequence),len(seqx),dot.x,dot.y,dot.length)
			pos, length = findLongestMatch( tables, sequence, seqx, self.ktup, self.window, self.mismatch, self.minmatch )

			#if pos==None:
				#print "No match found"
			#else:
				#print "pos=",pos,"length=",length
				
				#print seqx
				#print sequence[pos:pos+length]
			#print "="*40
			self.conservedstore[0].AddDot(dot.x,dot.y,length)
		
		end0=self.ProcEnd(0,None)
		end1=self.ProcEnd(1,None)
		translate = lambda seq: ''.join([{'A':'T','T':'A','G':'C','C':'G'}[a] for a in seq.upper()])
		for dotnum in range(len(backward)):
			print dotnum,"/",len(forward)
			dot = backward.GetDot(dotnum)
			seqx = self.seqstore[0][dot.x:dot.x+dot.length]
			#seqx = self.seqstore[0][end0-dot.x-dot.length:end0-dot.x]
			seqy = translate(self.seqstore[1][end1-dot.y-dot.length:end1-dot.y])[::-1]
			#seqy = self.seqstore[1][dot.y:dot.y+dot.length]
			#seqx = self.GetSubSequence(0,dot.x,dot.x+dot.length)
			#seqy = self.GetSubSequence(1,dot.y,dot.y+dot.length)
			#print seqx
			#print seqy
			assert seqx==seqy

			# now using the tables indexing, search for the longest match from this plot sequence.
			pos, length = findLongestMatch( tables, sequence, seqx, self.ktup, self.window, self.mismatch, self.minmatch )

			#if pos==None:
				#print "No match found"
			#else:
				#print "pos=",pos,"length=",length
				
				#print seqx
				#print sequence[pos:pos+length]
			#print "="*40
			self.conservedstore[1].AddDot(dot.x,dot.y,length)
			
		self.dotstore["conserved"]=self.conservedstore
		
		
	def CalculateDotStore(self, dimension=1, start=None, end=None, compstart=None, compend=None ):
		"""
		\brief performs the calculation of the dot plot for the entire sequences
		\details Calls the underlying libfreckle C code to calculate the dot plot for the entire sequence
		\warning can be slow with large sequences. 500000 vs 500000 can take a few hours to compute on 3 gigs of 32bit CPU
		\param dimension Which dimension to calculate the dotstore for (other dimension must be indexed)
		\param start the starting index number
		\param end the ending index number
		\todo make this method able to calculate sub dot plots.
		\todo implement ktuple size, window and mismatch
		\return (forward dotstore, reverse dotstore)
		"""
		assert(self.minmatch>=self.ktup)
		assert(self.ktup>=4)
		assert(self.window>=self.ktup)
		assert(dimension==0 or dimension==1)
		
		start=self.ProcStart(start)
		end=self.ProcEnd(dimension,end)
		compstart=self.ProcStart(compstart)
		compend=self.ProcEnd(1-dimension,compend)
		
		# make sure our other dimension is indexed in the specified region
		tables=self.tables[1-dimension][(compstart,compend)]
		assert(tables)		# other dimension should be indexed
		
		# assemble our comparison sequence
		compseq=decodeseq(self.GetSubSequence(dimension,start,end))
		
		# make a dotstore for this region
		dotstore=self.Compare(tables[3], tables[2], compseq, self.ktup, self.window, self.mismatch, self.minmatch)
		
		# and a reverse dotstore
		revdotstore=self.Compare(tables[3], tables[2], compseq[::-1], self.ktup, self.window, self.mismatch, self.minmatch)
		self.dotstore[ (dimension,start,end,compstart,compend) ] = (dotstore, revdotstore)
				
		# make sure the dotstore sizes are the same (and maximal)
		maxx=max(dotstore.GetMaxX(), revdotstore.GetMaxX())
		maxy=max(dotstore.GetMaxY(), revdotstore.GetMaxY())
		
		dotstore.SetMaxX(maxx)
		dotstore.SetMaxY(maxy)
		revdotstore.SetMaxX(maxx)
		revdotstore.SetMaxY(maxy)
		
		return (dotstore, revdotstore)
	
	def Compare(self,table,tableseq,compseq,ktup,window,mismatch,minmatch):
		return doComparison(table,tableseq,compseq,ktup,window,mismatch,minmatch)
	
	def Save(self,filename):
		"""
		\brief saves the dotplot to a file
		"""
		file=open(filename,"wb")
		
		# version numbers
		MAJOR=0
		MINOR=1
		
		# write the header to the file
		file.write(pack("4sii","_FDP",MAJOR,MINOR))
		
		# write the generating parameters
		# ktup, window, minmatch, mismatch
		file.write(pack("iiii",self.ktup, self.window, self.minmatch, self.mismatch))
		
		# write the generating filenames and filebounds.
		pickle.dump([self.filenames, self.sequencebounds, self.sequenceboundids,
			self.filebounds, self.size, self.globalfilebounds, self.globalsequencebounds], file)
		
		# write a datetime stamp
		file.write(pack("i",int(time())))
		
		# write number of dotstores
		file.write(pack("i",len(self.dotstore)))
		
		# write out each dot store
		for key in self.dotstore.keys():
			value=self.dotstore[key]
			
			#write key
			file.write(pack("iiiii",*key))
			
			# write forward and then backward store
			[ds.Save(file) for ds in value]
			
		# done. close the file
		file.close()
		

		
	def Load(self, filename):
		"""
		\brief loads the dotplot structure from a file
		"""
		file=open(filename,"rb")
		
		reader=lambda format: unpack(format, file.read( calcsize( format ) ) )
		
		#read the header
		head=reader("4s")[0]
		
		if head != "_FDP":
			file.close()
			raise DotPlotFileError, "Unknown file format"
		
		#read the version
		major,minor=reader("ii")
		
		#try to call the relevant reader
		try:
			result=eval("self.Load_%d_%d(file)"%(major,minor))
		except AttributeError, e:
			file.close()
			raise DotPlotFileError, "Cannot load version %d.%d of the freckle file format. Do you have the latest version?"%(major,minor)
		
		file.close()
		
		return result
	
	def Load_0_1(self,file):
		"""
		\brief load version 0.1 of the file format
		"""
		reader=lambda format: unpack(format, file.read( calcsize( format ) ) )
		
		# read the generating parameters
		# ktup, window, minmatch, mismatch
		self.ktup,self.window,self.minmatch,self.mismatch=reader("iiii")
		
		# write the generating filenames and filebounds.
		[self.filenames, self.sequencebounds, self.sequenceboundids,
			self.filebounds, self.size, self.globalfilebounds, self.globalsequencebounds]=pickle.load(file)
		
		# read the datetime stamp
		self.generatedon=reader("i")[0]
		
		# read number of dotstores
		numstores=reader("i")[0]
		
		# read each dotstore in
		self.dotstore={}
		
		for i in xrange(numstores):
			key=reader("iiiii")
			fds,rds=DotStore(),DotStore()
			fds.Load(file)
			rds.Load(file)
			self.dotstore[key]=(fds,rds)
		
		
	
	def IndexDotStores(self):
		"""
		\brief Indexes all the calculated DotStore
		"""
		# create indexes
		[[dots.CreateIndex() for dots in stores] for stores in self.dotstore.values()]
	
	def MakeAverageGrid(self,scale,x1=None,y1=None,x2=None,y2=None,storekey=None, conserved=False):
		"""
		\brief Calculates the reduced score grid
		\details Uses a fairly optimised algorithm to convert the calculated dot store into a grid of averaged values. These
		values represent the count of how many dots would be in that area.
		\todo allow sub windows to be averaged
		\param scale the scale value for the sizing of the grid. >1 to shrink. eg 10 means the final averaged grid
		will be 1/10th the size of the full dotplot.
		\param storekey the key to the dotstore. Of the form (dimension, start, end, compstart, compend)
		\return the DotGrid class
		"""
		if storekey==None:
			storekey=(1,0,self.GetSequenceLength(1),0,self.GetSequenceLength(0))
			
		key="conserved" if conserved else storekey
		
		assert(self.dotstore[key])
		self.scale=scale
		
		(dim,start,end,compstart,compend)=storekey
		
		width=compend-compstart
		height=end-start
		
		if x1==None and x2==None and y1==None and y2==None:
			x1=y1=0
			x2=width
			y2=height
				
		# forward and reverse grid
		grid=[DotGrid(),DotGrid()]
		
		[g.Calculate(dots,x1,y1,x2,y2,scale,self.window) for g,dots in zip(grid,self.dotstore[key])]
		grid[1].FlipInplace()
		grid[0].AddInplace(grid[1])
		self.grid[key]=grid[0]
		
		return grid[0]
		
	def MakeImage(self, storekey=None,major=None,minor=None,seqbound=(0,0,255),filebound=(255,0,0),alpha=24, conserved=False, invert=False):
		"""
		\brief Makes a DotPlot image from the averaged grid data
		"""
		if storekey==None:
			storekey=(1,0,self.GetSequenceLength(1),0,self.GetSequenceLength(0))
		
		key="conserved" if conserved else storekey
		
		string=self.grid[key].ToString()
		image=Image.fromstring("L", (self.grid[key].GetWidth(),self.grid[key].GetHeight()), string)
		
		if invert:
			image=ImageChops.invert(image)
		
		#image=image.convert("RGB")
		
		return image
		
	def CompileImage(self, storekey=None,major=None,minor=None,seqbound=(0,0,255),filebound=(255,0,0),alpha=24, dotimage=None, conservedimage=None):
		"""
		\brief Makes a DotPlot image from the averaged grid data
		"""
		if storekey==None:
			storekey=(1,0,self.GetSequenceLength(1),0,self.GetSequenceLength(0))
		
		if conservedimage:
			assert dotimage.size == conservedimage.size
		
		self.DrawBounds(dotimage,storekey[3],storekey[1],storekey[4],storekey[2], seqbound=tuple(list(seqbound)+[alpha]),filebound=tuple(list(seqbound)+[alpha]))
		image=self.AddAxis(dotimage,storekey[3],storekey[1],storekey[4],storekey[2],major=major,minor=minor,seqbound=seqbound,filebound=filebound)
		
		return image
	
	def AddAxis(self,image,x1,y1,x2,y2,major=None,minor=None, seqbound=(0,0,255), filebound=(255,0,0)):
		"""
		\brief add the image if an axis and sequence annotations to the dotplot image
		\detail takes the graph image and creates a new image with axis and annotations. the new image will be larger than the old
		\param image the pre calculated image of the dotplot
		\param x1 the left hand edge sequence offset
		\param y1 the top edge sequence offset
		\param x2 the right hand edge sequence offset
		\param y2 the bottom edge sequence offset
		\todo refactor this beheamoth
		\return the annotated image
		"""
		size=image.size
		font = getfont("/usr/share/fonts/truetype/freefont/FreeSans.ttf", 12)
		idfont = getfont("/usr/share/fonts/truetype/freefont/FreeSans.ttf", 10)
		titlefont = getfont("/usr/share/fonts/truetype/freefont/FreeSansBold.ttf", 14)
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
		
		# autocalculate our scale marks
		
		if major==None and minor!=None:
			major=10*minor
		elif minor==None and major!=None:
			minor=major/10
		elif minor==None and major==None:
			
			majors={
					0:	(10,1),
					50:	(500,50),
					100:	(10000,2000),
					500:	(50000,10000),
					5000:	(100000,10000),
					10000:	(100000,10000),
					100000:	(1000000,100000),
					1000000:	(10000000,1000000),
					10000000:	(100000000,10000000),
					100000000:	(1000000000,100000000),
					
				}
				
			
			location=sorted(majors.keys()+[self.scale]).index(self.scale)
			major,minor=majors[sorted(majors.keys())[location]]
		
		assert(major!=None)
		assert(minor!=None)
		
		MINORSCALE,MAJORSCALE=minor,major
		
		files=[[os.path.basename(a) for a in self.filenames[axis]] for axis in (0,1)]
		fileymaxlen=max([font.getsize(a)[0] for a in files[1]])
		filexmaxlen=max([font.getsize(a)[0] for a in files[0]])
		
		# find the side width
		ids=[]
		[[ids.append(a) for a in b] for b in self.sequenceboundids[0]]		# 0 is dimension (x)
		idxmaxlen=max([idfont.getsize(a)[0] for a in ids])
		
		# find the bottom height
		idys=[]
		[[idys.append(a) for a in b] for b in self.sequenceboundids[1]]		# 1 is dimension (y)
		idymaxlen=max([idfont.getsize(a)[0] for a in idys])
		
		bottomimagesize=(idxmaxlen+filexmaxlen+2*SIDEGUTTER, size[0])		
		
		# work out the largest sequence number
		lenseq1=x2-x1
		lenseq2=y2-y1
		digits=len(str(max(lenseq1,lenseq2)))			#this is the maximum number of digits in the sequence string
		
		axistextwidth,axistextheight=font.getsize("0"*digits)		#how wide the digit string is in pixels, and how high
		ticklength=10
		
		aisize=(size[0]+axistextwidth+ticklength+2+2*XOFF+idymaxlen+fileymaxlen+2*SIDEGUTTER,size[1]+2+2*YOFF+axistextheight+ticklength+titleheight+TITLEGAP+bottomimagesize[0]+YOFF)
		axisimage=Image.new("RGB",aisize,(255,255,255))
		
		dc = ImageDraw.Draw(axisimage,"RGBA")
		dc.rectangle( [(axistextwidth+ticklength+XOFF,axistextheight+ticklength+YOFF+titleheight+TITLEGAP),(aisize[0]-XOFF-1-idymaxlen-fileymaxlen-2*SIDEGUTTER,aisize[1]-YOFF-1-bottomimagesize[0]-YOFF)], fill=(0,0,0,255) )
		dc.rectangle( [(axistextwidth+ticklength+XOFF+2,axistextheight+ticklength+YOFF+titleheight+TITLEGAP+2),(aisize[0]-XOFF+1-idymaxlen-fileymaxlen-2*SIDEGUTTER,aisize[1]-YOFF+1-+bottomimagesize[0]-YOFF)], fill=(0,0,0,255) )
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
		seqys=[0]+self.seqybounds+[size[1]]
		seqyoff=1
		for file in boundids[axis]:
			for ids in file:
				idy=seqys[seqyoff]-(seqys[seqyoff]-seqys[seqyoff-1])/2
				
				textw,texth=idfont.getsize(ids)
				x=aisize[0]-XOFF-idymaxlen-fileymaxlen-SIDEGUTTER-SIDEGUTTER/2
				y=YOFF+titleheight+TITLEGAP+ticklength+axistextheight+idy-texth/2
				dc.text( (x,y), ids, fill=tuple(list(seqbound)+[255]), font=idfont)
				
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
			dc.text( (x,y), file, fill=tuple(list(filebound)+[255]), font=font)
			
			# draw the lines
			dc.line( [ (aisize[0]-XOFF-fileymaxlen-SIDEGUTTER, fbounds[fyoff-1]+tops+1), (aisize[0]-XOFF-fileymaxlen-2*SIDEGUTTER/3, fbounds[fyoff-1]+tops+1) ], fill=tuple(list(filebound)+[128]) )
			dc.line( [ (aisize[0]-XOFF-fileymaxlen-SIDEGUTTER, fbounds[fyoff]+tops-1), (aisize[0]-XOFF-fileymaxlen-2*SIDEGUTTER/3, fbounds[fyoff]+tops-1) ], fill=tuple(list(filebound)+[128]) )
			dc.line( [ (aisize[0]-XOFF-fileymaxlen-2*SIDEGUTTER/3, fbounds[fyoff-1]+tops+1), (aisize[0]-XOFF-fileymaxlen-2*SIDEGUTTER/3, fbounds[fyoff]+tops-1)], fill=tuple(list(filebound)+[128]))
			dc.line( [ (aisize[0]-XOFF-fileymaxlen-2*SIDEGUTTER/3,y+fnameh/2),(aisize[0]-XOFF-fileymaxlen-SIDEGUTTER/3,y+fnameh/2)],fill=tuple(list(filebound)+[128]))
			
			fyoff+=1
		
		################################################
		# horizontal annotations at the image bottom
		# we draw these in the regular way, then rotate them and paste them at the bottom
		################################################
		# find the side width
		ids=idys
		
		bottomimage=Image.new("RGB",bottomimagesize,(255,255,255))
		dc = ImageDraw.Draw(bottomimage,"RGBA")
		
		# size sequence annotations
		axis=0
		files=self.filenames
		boundids=self.sequenceboundids
		seqxs=[0]+self.seqxbounds+[size[0]]
		seqxoff=1
		for file in boundids[axis]:
			for ids in file:
				idx=seqxs[seqxoff]-(seqxs[seqxoff]-seqxs[seqxoff-1])/2
				
				textw,texth=idfont.getsize(ids)
				y=bottomimagesize[1]-idx-texth/2
				dc.text( (0,y), ids, fill=tuple(list(seqbound)+[255]), font=idfont)
				
				seqxoff+=1
		
		# file annotations
		axis=0
		files=[os.path.basename(a) for a in self.filenames[axis]]
		fbounds=[0]+[int(a/self.scale) for a in self.globalfilebounds[axis]]
		fyoff=1
		for file in files:
			idy=fbounds[fyoff]-(fbounds[fyoff]-fbounds[fyoff-1])/2
			
			fnamew,fnameh=font.getsize(file)
			x=bottomimagesize[0]-filexmaxlen
			tops=0
			y=bottomimagesize[1]-tops-idy-fnameh/2
			dc.text( (x,y), file, fill=tuple(list(filebound)+[255]), font=font)
			
			# draw the lines
			dc.line( [ (bottomimagesize[0]-YOFF-filexmaxlen-SIDEGUTTER, bottomimagesize[1]-(fbounds[fyoff-1]+tops+1)), (bottomimagesize[0]-YOFF-filexmaxlen-2*SIDEGUTTER/3, bottomimagesize[1]-(fbounds[fyoff-1]+tops+1)) ], fill=tuple(list(filebound)+[128]) )
			dc.line( [ (bottomimagesize[0]-YOFF-filexmaxlen-SIDEGUTTER, bottomimagesize[1]-(fbounds[fyoff]+tops-1)), (bottomimagesize[0]-YOFF-filexmaxlen-2*SIDEGUTTER/3, bottomimagesize[1]-(fbounds[fyoff]+tops-1)) ], fill=tuple(list(filebound)+[128]) )
			dc.line( [ (bottomimagesize[0]-YOFF-filexmaxlen-2*SIDEGUTTER/3, bottomimagesize[1]-(fbounds[fyoff-1]+tops+1)), (bottomimagesize[0]-YOFF-filexmaxlen-2*SIDEGUTTER/3, bottomimagesize[1]-(fbounds[fyoff]+tops-1))], fill=tuple(list(filebound)+[128]))
			dc.line( [ (bottomimagesize[0]-YOFF-filexmaxlen-2*SIDEGUTTER/3,(y+fnameh/2)),(bottomimagesize[0]-YOFF-filexmaxlen-SIDEGUTTER/3,y+fnameh/2)],fill=tuple(list(filebound)+[128]))
			
			fyoff+=1
		
		# rotate the image
		axisimage.paste(bottomimage.rotate(270),(axistextwidth+ticklength+2+XOFF,size[1]+2+2*YOFF+axistextheight+ticklength+titleheight+TITLEGAP))
			
		return axisimage
		
		
	
	def DrawBounds(self,image,xstart,ystart,xend,yend,filebound=(255,0,0,24),seqbound=(0,0,255,24)):
		"""
		draw the bounds as lines on the image
		"""
		#print "DrawBounds(",image,",",xstart,",",ystart,",",xend,",",yend,")"
		dc = ImageDraw.Draw(image,"RGBA")
		
		#dc.text((10, 25), "world", font=font,fill=(0,0,0,255))
		
		assert(xend>xstart)
		assert(yend>ystart)
		scale=self.scale
		
		#find sequence bounds and store their relative position to this image
		seqxbounds=[int(float(value-xstart)/scale) for value in reduce(lambda a,b:a+b,self.globalsequencebounds[0]) if value>=xstart and value<=xend]
		seqybounds=[int(float(value-ystart)/scale) for value in reduce(lambda a,b:a+b,self.globalsequencebounds[1]) if value>=ystart and value<=yend]
			
		self.seqybounds=seqybounds[:]
		self.seqxbounds=seqxbounds[:]
			
		#print "globalseqbounds[0]",len(self.globalsequencebounds[0][0])
		#print "globalseqbounds[1]",len(self.globalsequencebounds[1][0])
		#print "seqxbounds",len(self.seqxbounds)
		#print "seqybounds",len(self.seqybounds)
		
			
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
		
	def GetSubSequence(self, dimension=0, start=None, end=None):
		"""
		\brief return an assembled subsequence of the dotplot
		\detail return a sequence that is a combination of one or more sequences, from global seq offset 'start' to 'end'
		\param dimension the dimension to retrieve from. 0 for x. 1 for y
		\param start the start sequence number
		\param end the end sequence number
		\return a new sequence object
		"""
		assert(dimension==0 or dimension==1)		# atm only 2D
		
		# preprocess start and end
		start = self.ProcStart(start)
		end = self.ProcEnd(dimension,end)
		
		#print start,end,self.size,dimension
		
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



class LBDotPlot(DotPlot):
	"""Overrides the standard DotPlot class and replaces dot plot calculation with lbdot calculator"""
	
	def CreateTables(self, dimension=0, start=None, end=None):
		pass
	
	def CalculateDotStore(self, dimension=1, start=None, end=None, compstart=None, compend=None ):
		"""
		\brief performs the calculation of the dot plot for the entire sequences
		\details Calls the underlying libfreckle C code to calculate the dot plot for the entire sequence
		\warning can be slow with large sequences. 500000 vs 500000 can take a few hours to compute on 3 gigs of 32bit CPU
		\param dimension Which dimension to calculate the dotstore for (other dimension must be indexed)
		\param start the starting index number
		\param end the ending index number
		\todo make this method able to calculate sub dot plots.
		\todo implement ktuple size, window and mismatch
		\return (forward dotstore, reverse dotstore)
		"""
		assert(self.minmatch>=self.ktup)
		assert(self.ktup>=4)
		assert(self.window>=self.ktup)
		assert(dimension==0 or dimension==1)
		
		start=self.ProcStart(start)
		end=self.ProcEnd(dimension,end)
		compstart=self.ProcStart(compstart)
		compend=self.ProcEnd(1-dimension,compend)
		
		# assemble our comparison sequence
		compseq=decodeseq(self.GetSubSequence(dimension,start,end))
		tableseq=decodeseq(self.GetSubSequence(1-dimension,compstart,compend))
		
		# make a dotstore for this region
		dotstore,revdotstore=self.Compare(None, tableseq, compseq, self.ktup, self.window, self.mismatch, self.minmatch)
		self.dotstore[ (dimension,start,end,compstart,compend) ] = (dotstore, revdotstore)
				
		# make sure the dotstore sizes are the same (and maximal)
		maxx=max(dotstore.GetMaxX(), revdotstore.GetMaxX())
		maxy=max(dotstore.GetMaxY(), revdotstore.GetMaxY())
		
		#print maxx, maxy
		
		dotstore.SetMaxX(maxx)
		dotstore.SetMaxY(maxy)
		revdotstore.SetMaxX(maxx)
		revdotstore.SetMaxY(maxy)
		
		
		return (dotstore, revdotstore)
	
	def Compare(self,table,tableseq,compseq,ktup,window,mismatch,minmatch):
		return doFastComparison(tableseq,compseq,ktup,window,mismatch,minmatch)


#	
# test suite
#
import unittest, random

class TestDotPlot(unittest.TestCase):
	NUMAXIS=2
	
	def setUp(self):
		from glob import glob
		self.filelist=glob("TestFastaFiles/*.fasta")
	
	def testBounds(self):
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
				
	def testGetSubSequence(self):
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
		
	def testCompareSequences(self):
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
					
			
	def testMakeDotPlotSubImage(self):
		dp=DotPlot(self.filelist, self.filelist)
		
		fullseq=dp.AssembleFullSequence(0)
		seqlen=len(fullseq)
		
		image=dp.MakeDotPlotSubImage(fullseq,fullseq)
		image.save("test.png")
		
	def testBigPlot(self):
		dp=DotPlot(self.filelist, self.filelist)
		
		image=dp.MakeDotPlotSubImage(dp.GetSubSequence(0,0,1024),dp.GetSubSequence(1,0,512),-5)
		image.save("test.png")
		
	def testPlotWithBorders(self):
		dp=DotPlot(self.filelist, self.filelist)
		
		image=dp.DrawDotPlot(0,0,window=4)
		image.save("test.png")
		
	def testCreateTables(self):
		dp=DotPlot(self.filelist, self.filelist)
		
		print "tables"
		tables=dp.CreateTables()
		
		print "dotstore"
		dotstore=dp.CalculateDotStore()
		
		print "index"
		dp.IndexDotStores()
		
if __name__ == '__main__':
    unittest.main()
