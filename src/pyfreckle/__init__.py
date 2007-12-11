#
# Module to use libfreckle in python
#

"""\
pyfreckle is the python bindings of libfreckle, to aid in the computation and construction of dot plots
"""

from ctypes import cdll
import os.path

# dynamically locate our .so file in the same directory as this module
__library_c_source__=os.path.join(os.path.dirname(__file__),"libfreckle.so")
lib=cdll.LoadLibrary(__library_c_source__)

# import the modules into this namespace
from DotGrid import DotGrid
from DotStore import DotStore

# set a static class variable that is the library
DotGrid.lib=lib
DotStore.lib=lib

# now our base library functions
def buildMappingTables( sequence, ktuplesize ):
	return lib.buildMappingTables(sequence, ktuplesize)

def makeDotComparison(seq1, seq2, ktuplesize=4, window=10, mismatch=0, minmatch=4):
	return DotStore(lib.makeDotComparison(seq1,seq2,ktuplesize,window,mismatch,minmatch))
