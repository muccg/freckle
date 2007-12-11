#!/usr/bin/env python

#
# PyDotter
# ========
# PyDotter is a dot plot generator. It is designeed to plot multiple sequences against multiple
# sequences with the final graph having the originating sequences identified.

import sys, getopt
from Bio import SeqIO
from numpy import array, zeros, int32
from PIL import Image

# switch this to True to get debug messages
DEBUG=True

def usage():
	"""print the usage of the programme"""
	print "Usage:"
	print " %s [-h] -x xseq1 [-x xseq2 -x xseq3 ...] -y yseq1 [-y yseq2 -y yseq3 ...] -o outputfile [-s imagesize]"%sys.argv[0]
	print "-h\tthis helpful blurb"
	print "-x\ta FASTA file for a sequence to be displayed on the x axis"
	print "-y\ta FASTA file for a sequence to be displayed on the y axis"
	print "-o\tthe image file to write the output to. Use '-' for stdout"
	print "-s\tthe output image size. Should be expressed as integers seperated by an 'x'. eg. '800x600'"
	print
	print "Examples:"
	print " %s -x seq1.fasta -x seq2.fasta -x seq3.fasta -y seq4.fasta -y seq5.fasta -o graph.png -s 1024x768"%sys.argv[0]
	
def parseopts():
	"""parse the command line options and return them"""
	
	# default parameters
	xseq=[]
	yseq=[]
	outfile="output.png"
	imagesize=512
	
	#our getopt definition strings
	shortopts="hx:y:o:s:"
	longopts=["help","xfile=","yfile=","output=","size="]
	
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
			
		elif o in ("-s","--size"):
			# try and parse the size string
			try:
				size=parsesizestring(a)
			except Exception, e:
				print "ERROR: cannot parse size string!",str(e)
				usage()
				sys.exit(3)
				
	return xseq, yseq, outfile, imagesize
	
class ParseException(Exception):
	"""An error parsing the size string"""
	
def parsesizestring(string):
	# extract any unallowed characters with a list comp.
	unallowed=[a for a in string.lower() if a not in "0123456789x"]
	if len(unallowed):
		raise ParseException, "Unallowed characters in string"
	
	#split on the 'x'. There should only be two numbers
	numbers=string.lower().split('x')
	if len(numbers)!=2:
		raise ParseException, "Too many parameters in string"
	
	return tuple([int(x) for x in numbers])
	
def main():
	xseqfiles,yseqfiles,outfile,imagesize=parseopts()
	
	if DEBUG:
		print "xsequences:",xseqfiles
		print "ysequences:",yseqfiles
		print "outfile:",outfile
		print "imagesize",imagesize
	
	from DotPlot import DotPlot
	
	plot=DotPlot(xseqfiles,yseqfiles)
	
	print "Size = %d x %d"%(plot.GetSequenceLength(0),plot.GetSequenceLength(1))
	
	## try to conserve memory.
	## save the boundaries between files and between sequences.
	#xseqbounds,yseqbounds=[[[len(x.seq) for x in SeqIO.parse(open(file),"fasta")] for file in seqfiles] for seqfiles in xseqfiles,yseqfiles]

	## save the beginning offsets of each file
	#xfileoff,yfileoff=[[sum(c) for c in bounds] for bounds in [xseqbounds,yseqbounds] ]
	
	## how big is our full matrix
	#xsize,ysize=[sum([a for a in fileoff]) for fileoff in xfileoff,yfileoff]
	
	#print xseqbounds,yseqbounds
	#print "---"
	#print xfileoff,yfileoff
	#print "---"
	#print xsize,ysize
	
	#sys.exit()
	
	

	
	
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