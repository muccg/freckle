#
# Module to use libfreckle in python
#

"""\
pyfreckle is the python bindings of libfreckle, to aid in the computation and construction of dot plots
"""

from ctypes import *
import os

# dynamically locate our .so file in the same directory as this module
try:
	import sys
	if sys.frozen:
		# we are frozen in a package. Grab it from our unpacked directory
		__library_c_source__=os.path.join(os.environ['LD_LIBRARY_PATH'].split(':')[0],"libfreckle.so")
		lib=cdll.LoadLibrary(__library_c_source__)
except AttributeError, e:
	OPATHS=['./','/usr/local/lib','/usr/lib','/lib',os.path.dirname(__file__)]
	PATH=OPATHS[::-1]
	lib=None
	while lib==None:
		if len(PATH)==0:
			#no library could be found
			raise Exception, "libfreckle.so could not be found. tried " + ','.join(OPATHS)
		__library_c_source__=os.path.join(PATH.pop(),"libfreckle.so")
		try:
			lib=cdll.LoadLibrary(__library_c_source__)
		except OSError, e:
			pass
		

# import the modules into this namespace
from DotGrid import DotGrid
from DotStore import DotStore

# set a static class variable that is the library
DotGrid.lib=lib
DotStore.lib=lib

# set vairables
lib.Bases=c_char_p.in_dll(lib, "Bases")
lib.Aminos=c_char_p.in_dll(lib, "Aminos")

class c_void(Structure):
    # c_void_p is a buggy return type, converting to int, so
    # POINTER(None) == c_void_p is actually written as
    # POINTER(c_void), so it can be treated as a real pointer.
    _fields_ = [('dummy', c_int)]
    
# set passing and return types where needed
lib.buildMappingTables.argtypes = [POINTER(c_char), c_int, POINTER(c_char)]
lib.buildMappingTables.restype = POINTER(c_void)
lib.doComparison.argtypes=[POINTER(c_void), POINTER(c_char), POINTER(c_char), c_int, c_int, c_int, c_int, POINTER(c_char)]
lib.DotStoreToBuffer.restype = POINTER(c_int)
lib.DotStoreFromBuffer.argtypes = [ POINTER(c_void), POINTER(c_int) ]
lib.NewDotStore.restype=POINTER(c_void)
lib.DotStoreGetMaxX.restype=c_int
lib.DotStoreGetMaxY.restype=c_int
lib.DotStoreSetMaxX.argtypes=[c_void_p,c_int]
lib.DotStoreSetMaxY.argtypes=[c_void_p,c_int]
lib.DotGridToString.argtypes=[POINTER(c_void)]
lib.DotGridToString.restype=POINTER(c_void)
lib.NewDotGrid.argtypes=[]
lib.NewDotGrid.restype=POINTER(c_void)

class c_pointers(Structure):
	_fields_ = [ ('forward', POINTER(c_void)),('reverse',POINTER(c_void))]

lib.DoFastComparison.restype=POINTER(c_pointers)

# now our base library functions
#def buildMappingTables( sequence, ktuplesize ):
	#return lib.buildMappingTables(sequence, ktuplesize)

def makeDotComparison(seq1, seq2, ktuplesize=4, window=10, mismatch=0, minmatch=4):
	return DotStore(lib.makeDotComparison(seq1,seq2,ktuplesize,window,mismatch,minmatch))

def buildMappingTables( sequence, ktuplesize, alphabet=lib.Bases ):
	return lib.buildMappingTables(sequence, ktuplesize, alphabet)

#DotStore *doComparison(int **tables, const char *tablesequence, const char *newsequence, int ktuplesize, int window, int mismatch, int minmatch, const char *bases=Bases );
def doComparison(tables, tabseq, newseq, ktup, window, mismatch, minmatch, bases=lib.Bases):
	return DotStore(lib.doComparison(tables,tabseq,newseq,ktup,window,mismatch,minmatch,bases))
	
def doFastComparison(seq1, seq2, ktuplesize=4, window=10, mismatch=0, minmatch=4):
	results=lib.DoFastComparison(seq1,seq2,len(seq1),len(seq2),window,mismatch,0,ktuplesize)
	return DotStore(results.contents.forward),DotStore(results.contents.reverse)

	
