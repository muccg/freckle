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
		






import unittest

class TestSequenceFunctions(unittest.TestCase):
	def setUp(self):
		self.seq1 = "GCGGGTACTGATATACTCATGATTATACCGCGCGGTTGTGTGAATTAATATCAACACCACAAAAGAGAGGAGGACTTCCTCTCTCTCTCTAACACCAATATATCCGGCCGGTTG"
		self.seq2 = "ATCGACGTATAGATTTTTCCACAGCGCCAAACTCTTCTATCACTCATGACTGACTGTGTCATGACTGATTATATATATCTCTCTTCTCATATATCATACT"
		self.freckle=Freckle()
	
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
			
	def testAccuracy(self):
		self.freckle.Compare(self.seq1, self.seq2,ktuple=2,window=4,minmatch=4,mismatch=0)		#accuracte matching of length at least 4
		length=len(self.freckle)
		for i in xrange(length):
			x,y,length=self.freckle[i]
			print i,x,y,length
			print self.seq1[x:x+length],"<->",self.seq2[y:y+length]
			self.assert_(self.seq1[x:x+length]==self.seq2[y:y+length])
	
	
if __name__ == '__main__':
	unittest.main()
	
	
