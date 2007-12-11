from ctypes import *

class Dot(Structure):
	"""Wrapper for C struct Dot.
	
	struct
	{
		int x;
		int y;
		int length;
	};
	"""
	_fields_ = [ 	("x", c_int),
			("y", c_int),
			("length", c_int)	]
	
	
	def __repr__(self):
		old=Structure.__repr__(self)				#<pyfreckle.Dot.Dot object at 0xb7dc7104>
		parts=old[1:-1].split(" ")
		return "<%s (%d,%d,%d) at %s>"%(parts[0],self.x,self.y,self.length,parts[-1])

