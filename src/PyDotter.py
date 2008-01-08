#!/usr/bin/env python

#
# PyDotter
# ========
# PyDotter is a dot plot generator. It is designeed to plot multiple sequences against multiple
# sequences with the final graph having the originating sequences identified.

import sys, getopt
from Bio import SeqIO
from numpy import array, zeros, int32
from PIL import Image, ImageDraw
from time import time

# switch this to True to get debug messages
DEBUG=False

def usage():
	"""print the usage of the programme"""
	print "Usage:"
	print " %s [-h] -x xseq1 [-x xseq2 -x xseq3 ...] -y yseq1 [-y yseq2 -y yseq3 ...] -o outputfile [-s imagesize]"%sys.argv[0]
	print "-h\tthis helpful blurb"
	print "-x\ta FASTA file for a sequence to be displayed on the x axis"
	print "-y\ta FASTA file for a sequence to be displayed on the y axis"
	print "-o\tthe image file to write the output to."
	print "-s\tthe output image size. Should be expressed as integer. Determines the longest side of the grid. real image will be slightly larger."
	print "-k\tktuple size for tokenisation"
	print "-w\twindow size for comparing mismatches"
	print "-m\tminimum match length to be included"
	print "-d\tnumber of mismatches allowed per window of sequence"
	print "-S\tsave the computed dotplot as the specified file"
	print "-L\tload the precomputed dotplot from the specified file"
	print "-M\toverride the automatic major tick seperation with this value"
	print "-T\toverride the automatic minor tick seperation with this value"
	print
	print "Examples:"
	print " %s -x seq1.fasta -x seq2.fasta -x seq3.fasta -y seq4.fasta -y seq5.fasta -o graph.png -s 1024"%sys.argv[0]
	
def parseopts():
	"""parse the command line options and return them"""
	
	# default parameters
	xseq=[]
	yseq=[]
	outfile="output.png"
	imagesize=512
	ktup=8
	window=8
	minmatch=8
	mismatch=0
	savefile=None
	loadfile=None
	major=None
	minor=None
	
	#our getopt definition strings
	shortopts="hx:y:o:s:k:w:m:d:S:L:M:T:"
	longopts=["help","xfile=","yfile=","output=","size=","ktup=","window=","minmatch=","mismatch=","save=","load=","major=","minor="]
	
	try:
		opts,args=getopt.getopt(sys.argv[1:],shortopts,longopts)
	except getopt.GetoptError:
		usage()
		sys.exit(2)
		
	if len(args):
		print "ERROR: command line has extrenuous arguments"
		usage()
		sys.exit(1)
		
	for o,a in opts:
		if o in ("-h", "--help"):
			usage()
			sys.exit()
			
		elif o in ("-x","--xfile"):
			xseq.append(a)
			
		elif o in ("-y","--yfile"):
			yseq.append(a)
			
		elif o in ("-o","--output"):
			outfile=a
			
		elif o in ("-k","--ktup"):
			ktup=int(a)
			
		elif o in ("-w","--window"):
			window=int(a)
			
		elif o in ("-m","--minmatch"):
			minmatch=int(a)
			
		elif o in ("-d","--mismatch"):
			mismatch=int(a)
			
		elif o in ("-S","--save"):
			savefile=str(a)
			
		elif o in ("-L","--load"):
			loadfile=str(a)	
			
		elif o in ("-M","--major"):
			major=int(a)
			
		elif o in ("-T","--minor"):
			minor=int(a)
			
		elif o in ("-s","--size"):
			# try and parse the size string
			try:
				imagesize=int(a)
			except Exception, e:
				print "ERROR: cannot parse size string!",str(e)
				usage()
				sys.exit(3)
				
	return xseq, yseq, outfile, imagesize,ktup,window,minmatch,mismatch,savefile,loadfile,major,minor
	
def main():
	xseqfiles,yseqfiles,outfile,imagesize,ktup,window,minmatch,mismatch,savefile,loadfile,major,minor=parseopts()
	
	if DEBUG:
		print "xsequences:",xseqfiles
		print "ysequences:",yseqfiles
		print "outfile:",outfile
		print "imagesize",imagesize
		print "ktuple",ktup
		print "window",window
		print "minmatch",minmatch
		print "mismatch",mismatch
		print "loadfile",loadfile
		print "savefile",savefile
		print "major",major
		print "minor",minor
	
	from DotPlot import DotPlot
	
	plot=DotPlot(xseqfiles,yseqfiles,ktup, window, minmatch, mismatch)
	
	if loadfile!=None:
		#load the dotstore from a previous run
		print "Loading dotplot",loadfile,"..."
		t=time()
		plot.Load(loadfile)
		print "done in",time()-t,"seconds"
	else:
		#calculate it
		print "create tables"
		t=time()
		plot.CreateTables()
		print time()-t,"seconds"
	
		print "calculating dotplot"
		t=time()
		plot.CalculateDotStore()
		print time()-t,"seconds"
	
	
	if loadfile==None and savefile!=None:
		#save the dotstore
		print "Saving dotplot"
		plot.Save(savefile)
		
	
	xsize,ysize=plot.GetSequenceLength(0),plot.GetSequenceLength(1)
	
	# work out scale and thus final image width and height
	longest=(xsize>ysize) and xsize or ysize
	
	from math import ceil
	
	scale=ceil(float(longest)/float(imagesize))
	scale=float(longest)/float(imagesize)
	
	xoutput=float(xsize)/scale
	youtput=float(ysize)/scale
	
	print "indexing"
	t=time()
	plot.IndexDotStores()
	print time()-t,"seconds"
	
	print "averaging"
	t=time()
	plot.MakeAverageGrid(scale)
	print time()-t,"seconds"
	
	print "making image"
	t=time()
	image=plot.MakeImage(major=major,minor=minor)
	print time()-t,"seconds"
	
	image.save(outfile)
	
	sys.exit(0)
	
	

	
	
def dummy():
	
	# save the 
	
	if DEBUG:
		print "Matrix size: %d x %d"%(xsize,ysize)
	
	windowsize=1024*16
	window=zeros( (windowsize,windowsize), bool )
	
	destination=Image.new("RGB", imagesize, (255,255,255))
	
	for winx in xrange((xsize/windowsize)+1):
		for winy in xrange((ysize/windowsize)+1):
			if DEBUG:
				print "processing window %d,%d"%(winx,winy)
				
				seqx=getsubseq(xseqfiles,winx*windowsize)
				print seqx
				
				
def getsubseq(filelist,start,end):
	pass
				
		
	# create our whole array
	#points=zeros( (xsize>>2,ysize), int32 )
	
	
	
	#handle = open(xseqfiles[0])
	#for seq_record in SeqIO.parse(handle, "fasta") :
		#print seq_record.id
		#print seq_record.seq
		#print len(seq_record.seq)
	#handle.close()
		

if __name__ == "__main__":
    main()