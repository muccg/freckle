#
# tests the ctypes python importation of libfreckle
#

from ctypes import *
lib=CDLL("./libfreckle.so")
lib.getInfo(None)
