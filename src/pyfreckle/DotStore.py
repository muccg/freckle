
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
	
	def SetMaxX(self,x):
		self.lib.DotStoreSetMaxX(self.dotstore, x)
		
	def SetMaxY(self,y):
		self.lib.DotStoreSetMaxY(self.dotstore, y)
		
	def GetMaxX(self):
		return self.lib.DotStoreGetMaxX(self.dotstore)
		
	def GetMaxY(self):
		return self.lib.DotStoreGetMaxY(self.dotstore)
		
	
	
	
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
	
	def CreateIndex(self):
		self.lib.DotStoreCreateIndex(self.dotstore)
		
	def DestroyIndex(self):
		self.lib.DotStoreDestroyIndex(self.dotstore)
		
	def ToBuffer(self):
		return self.lib.DotStoreToBuffer(self.dotstore)
	
	def FromBuffer(self, buffer):
		self.lib.DotStoreFromBuffer(self.dotstore, buffer)
	
	def FreeBuffer(self, buffer):
		self.lib.FreeIntBuffer(buffer)
		
	def Filter(self, minmatch):
		return DotStore(self.lib.DotStoreFilter(self.dotstore, minmatch))
	
	def ToString(self):
		import struct
		
		# get our buffer and length
		buff=self.ToBuffer()
		length=self.lib.DotStoreBufferSize(self.dotstore,buff)
		
		# write it to the stream
		data=struct.pack("%di"%(length), *buff[:length])
		
		# free the buffer
		self.FreeBuffer(buff)
		
		return data
	
	def FromString(self, string):
		import struct
		
		data=string[:struct.calcsize("3i")]
		maxx,maxy,length=struct.unpack("3i",data)
		
		# turn into a shitload of ints. We don't need to free this because it was created in python
		array=apply( c_int*(length+3), struct.unpack("%di"%(length+3),string) )
		
		print len(array)
		
		#sys.exit(0)
		
		self.FromBuffer(array)
		
	
	def Load(self,stream):
		import struct
		
		# read the first three elements from the stream
		data=stream.read(struct.calcsize("3i"))
		maxx,maxy,length=struct.unpack("3i",data)
		
		#print maxx,maxy,length,struct.calcsize("%di"%length)
		
		# read length more
		data+=stream.read(struct.calcsize("%di"%length*3))
		
		#print "read",len(data)
		
		#array=c_int*(length+3)
		#print array
		
		# turn into a shitload of ints. We don't need to free this because it was created in python
		array=apply( c_int*(length*3+3), struct.unpack("%di"%(length*3+3),data) )
		
		#print "applied",len(array)
		
		self.FromBuffer(array)
	
	def Save(self,stream):
		stream.write(self.ToString())
	
	
	