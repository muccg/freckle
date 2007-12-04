#
# A thin wrapper to make accessing the C++ code in libfreckle a whole lot easier
# Also comes with a test harness!
#

LIBPATH="./libfreckle.so"

from ctypes import *

class Freckle:
	def __init__(self):
		self.lib=CDLL(LIBPATH)
		self.dotstore=None
		
	def __del__(self):
		if self.dotstore!=None:
			self.lib.FreeDotStore(self.dotstore)
		
	def Compare(self, seq1, seq2, ktuple=4, window=10, mismatch=1, minmatch=5):
		self.dotstore=self.lib.makeDotComparison(seq1, seq2, ktuple, window, mismatch, minmatch)
		
	def __getitem__(self, item):
		assert(self.dotstore!=None)
		lib=self.lib
		ds=self.dotstore
		num=lib.GetNumDots(ds)
		if item>=num or item< -num:
			raise IndexError, "Index out of range"
		if item>=0:
			return (lib.GetDotX(ds,item), lib.GetDotY(ds,item), lib.GetDotLength(ds,item))
		else:
			return self.__getitem__(num+item)		#recurse as a positive
		
	def __len__(self):
		assert(self.dotstore!=None)
		return self.lib.GetNumDots(self.dotstore)
	
	def MakeList(self):
		"""Turn the results into a big python list
		WARNING: could be a big mistake with a big results set
		"""
		return [self[i] for i in xrange(len(self))]
	
	def MakeGenerator(self):
		"""Turn the results into a generator"""
		return (self[i] for i in xrange(len(self)))
		






import unittest
import random

class TestSequenceFunctions(unittest.TestCase):
	def setUp(self):
		self.seq1 = "GCGGGTACTGATATACTCATGATTATACCGCGCGGTTGTGTGAATTAATATCAACACCACAAAAGAGAGGAGGACTTCCTCTCTCTCTCTAACACCAATATATCCGGCCGGTTG"
		self.seq2 = "ATCGACGTATAGATTTTTCCACAGCGCCAAACTCTTCTATCACTCATGACTGACTGTGTCATGACTGATTATATATATCTCTCTTCTCATATATCATACT"
		self.freckle=Freckle()
		
	def tearDown(self):
		del self.freckle
		self.freckle=None
	
	def testCompare(self):
		self.freckle.Compare(self.seq1, self.seq2)
		self.assert_(self.freckle.dotstore)
		
	def testLength(self):
		self.freckle.Compare(self.seq1, self.seq2)
		self.assert_(len(self.freckle))
		
	def testGetItem(self):
		self.freckle.Compare(self.seq1, self.seq2)
		length=len(self.freckle)
		for i in xrange(length):
			# negative indexes work
			self.assert_(self.freckle[i]==self.freckle[i-length])
			
			# the type is a 3-tuple
			self.assert_(len(self.freckle[i])==3)
			self.assert_(type(self.freckle[i])==tuple)
			
	def testPrecision(self):
		"""Makes sure that each result is a true result"""
		self.freckle.Compare(self.seq1, self.seq2,ktuple=2,window=4,minmatch=4,mismatch=0)		#accuracte matching of length at least 4
		length=len(self.freckle)
		for i in xrange(length):
			x,y,length=self.freckle[i]
			#print i,x,y,length
			#print self.seq1[x:x+length],"<->",self.seq2[y:y+length]
			self.assert_(self.seq1[x:x+length]==self.seq2[y:y+length])
			
	def testAccuracy(self):	
		"""Makes sure that EVERY true result is represented"""
		seq1=self.seq1
		seq2=self.seq2
		output=[]
		for s1 in xrange(len(seq1)):
			for s2 in xrange(len(seq2)):
				matchlen=0
				try:
					while seq1[s1+matchlen]==seq2[s2+matchlen]:
						matchlen+=1
				except IndexError, e:
					pass
				if matchlen:
					output.append( (s1,s2,matchlen) )
					
		# only those of match length at least 4
		outfiltered=[out for out in output if out[2]>=4]
		self.freckle.Compare(seq1, seq2, ktuple=2,window=4,minmatch=4,mismatch=0)
		freckleout=self.freckle.MakeList()
		
		for item in outfiltered:
			for freckitem in freckleout:
				#print "<",item,",",freckitem,">"
				if item[:]==freckitem[:]:
					#print "removing",freckitem
					freckleout.remove(freckitem)
		
		self.assert_( len(freckleout)==0 )
		
	def testBigPrecision(self):
		"""Make two enormous strings and compute. Then make sure each result is correct"""
		bases="ACGT"
		seq1="".join([random.choice(bases) for i in xrange(random.randint(10000,20000))])
		seq2="".join([random.choice(bases) for i in xrange(random.randint(10000,20000))])
		
		self.freckle.Compare(seq1, seq2, ktuple=4,window=4,minmatch=4,mismatch=0)
		length=len(self.freckle)
		matches,misses=0,0
		for i in xrange(length):
			x,y,length=self.freckle[i]
			self.assert_(seq1[x:x+length]==seq2[y:y+length])
			
	def testTime(self):
		bases="ACGT"
		seq1="".join([random.choice(bases) for i in xrange(350000)])
		seq2="".join([random.choice(bases) for i in xrange(300000)])
		
		import time
		t1=time.time()
		self.freckle.Compare(seq1, seq2, ktuple=12,window=12,minmatch=4,mismatch=0)
		print "k=12, w=12, min=4, mis=0 =>",time.time()-t1
		
		t1=time.time()
		self.freckle.Compare(seq1, seq2, ktuple=8,window=8,minmatch=8,mismatch=0)
		print "k=8, w=8, min=8, mis=0 =>",time.time()-t1
		
		t1=time.time()
		self.freckle.Compare(seq1, seq2, ktuple=4,window=4,minmatch=4,mismatch=0)
		print "k=4, w=4, min=4, mis=0 =>",time.time()-t1
		
		t1=time.time()
		self.freckle.Compare(seq1, seq2, ktuple=2,window=2,minmatch=2,mismatch=0)
		print "k=2, w=2, min=2, mis=0 =>",time.time()-t1
		
		
		
		
	
	
if __name__ == '__main__':
	unittest.main()
	
	
