class DotGrid:
	def __init__(self, dotgrid=None):
		if dotgrid==None:
			self.dotgrid=self.lib.NewDotGrid()
		else:
			self.dotgrid=dotgrid
		
	def __del__(self):
		assert(self.dotgrid)
		self.lib.DelDotGrid(self.dotgrid)
		self.dotgrid=None
		
	def Create(self,x,y):
		self.lib.DotGridCreate(self.dotgrid,x,y)
		
	def GetWidth(self):
		return self.lib.DotGridWidth(self.dotgrid)
	
	def GetHeight(self):
		return self.lib.DotGridHeight(self.dotgrid)
	
	def GetSize(self):
		return self.lib.DotGridGetSize(self.dotgrid)
	
	def GetPoint(self,x,y):
		return self.lib.DotGridGetPoint(self.dotgrid,x,y)
	
	def GetData(self,pos):
		return self.lib.DotGridGetData(self.dotgrid,pos)
	
	def SetPoint(self,x,y,data):
		self.lib.DotGridSetPoint(self.dotgrid,x,y,data)
	
	def GetMax(self):
		return self.lib.DotGridGetMax(self.dotgrid)
	
	def GetMin(self):
		return self.lib.DotGridGetMin(self.dotgrid)
	
	def ToString(self):
		string=self.lib.DotGridToString(self.dotgrid)
		return string_at(string,self.GetSize())
		
	def Calculate(self, source, x1, y1, x2, y2, scale, window):
		self.lib.DotGridCalculate(self.dotgrid, source, x1, y1, x2, y2, scale, window)
		
	def AddInplace(self, dotgrid):
		self.lib.DotGridAddInplace(self.dotgrid, dotgrid.dotgrid)
		
	def FlipInplace(self):
		self.lib.DotGridFlipInplace(self.dotgrid)





import unittest
import random

class TestSequenceFunctions(unittest.TestCase):
	def setUp(self):
		self.dg=DotGrid()
		
	def tearDown(self):
		del self.dg
		self.dg=None
	
	def testPass(self):
		self.assert_(True)
		
	def testCreate(self):
		self.dg.Create(100,200)
		self.assert_(self.dg.GetWidth()==100)
		self.assert_(self.dg.GetHeight()==200)
		self.assert_(self.dg.GetSize()==20000)
		
	def testGetSetPoint(self):
		self.dg.Create(10,20)
		for y in xrange(20):
			for x in xrange(10):
				self.dg.SetPoint(x,y,x+y)
				self.assert_(self.dg.GetPoint(x,y)==x+y)
		
	def testGetMaxMin(self):
		self.dg.Create(10,20)
		for y in xrange(20):
			for x in xrange(10):
				self.dg.SetPoint(x,y,x+y)
		
		self.assert_(self.dg.GetMax()==19+9)
		self.assert_(self.dg.GetMin()==0)
		
	def testToString(self):
		self.dg.Create(16,16)
		for y in xrange(16):
			for x in xrange(16):
				self.dg.SetPoint(x,y,x*y)
				
		self.assert_(self.dg.GetMax()==225)
		self.assert_(self.dg.GetMin()==0)
		
		string=self.dg.ToString()
		self.assert_(len(string)==256)
		
		# should be scaled between 0 and 255
		self.assert_(ord(string[0])==0)
		self.assert_(ord(string[-1])==255)
		
	def testAddInplace(self):
		self.dg.Create(16,16)
		for y in xrange(16):
			for x in xrange(16):
				self.dg.SetPoint(x,y,x*y)
				
		dg2=DotGrid()
		dg2.Create(16,16)
		for y in xrange(16):
			for x in xrange(16):
				dg2.SetPoint(x,y,x+y)
				
		self.dg.AddInplace(dg2)
		
		for y in xrange(16):
			for x in xrange(16):
				self.assert_(self.dg.GetPoint(x,y)==x+y+x*y)
				
		del dg2
		
	def testFlipInplace(self):
		self.dg.Create(16,16)
		for y in xrange(16):
			for x in xrange(16):
				self.dg.SetPoint(x,y,x*y)
				
		self.dg.FlipInplace()
		
		for y in xrange(16):
			for x in xrange(16):
				self.assert_(self.dg.GetPoint(x,y)==x*(15-y))
				
	
		
	
	
if __name__ == '__main__':
	from ctypes import *
	DotGrid.lib=cdll.LoadLibrary("./libfreckle.so")
	unittest.main()
	
	

