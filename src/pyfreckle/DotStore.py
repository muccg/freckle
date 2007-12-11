
import sys
from ctypes import *

from Dot import Dot

class DotStore:
	def __init__(self, dotstore=None):
		if dotstore==None:
			self.dotstore=self.lib.NewDotStore()
		else:
			self.dotstore=dotstore
			
		self.lib.DotStoreGetDot.restype=POINTER(Dot)
		
	def __del__(self):
		assert(self.dotstore)
		self.lib.DelDotStore(self.dotstore)
	
	def __len__(self):
		return self.lib.DotStoreGetNumDots(self.dotstore)
	
	def __getitem__(self, item):
		if type(item)==slice:
			start,stop,step=item.start,item.stop,item.step
			return [self.GetDot(i) for i in range(start, ((stop==sys.maxint) and len(self) or stop), ((step==None)and 1 or step))] 
		else:
			return self.GetDot(item)
	
	def GetDot(self, index):
		if index<0:
			index=len(self)+index
		assert(index>=0 and index<len(self))
		return self.lib.DotStoreGetDot(self.dotstore,index).contents
	
	def GetDotX(self, index):
		return self.lib.DotStoreGetDotX(self.dotstore,index)
		
	def GetDotY(self,index):
		return self.lib.DotStoreGetDotY(self.dotstore,index)
	
	def GetDotLength(self,index):
		return self.lib.DotStoreGetDotLength(self.dotstore,index)
	
	
	
	