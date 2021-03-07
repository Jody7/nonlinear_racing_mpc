import math, time
import numpy as np
from bisect import bisect_right

class PiecePoly():
	def horner(self, lst, x):
		total = 0
		for a in lst:
			total = total*x+a
		return total
	def __init__(self):
		form = 'pp'
		breaks = []
		coefs = []
		pieces = 0
		order = 0
		dim = 0
	def eval(self, x):
		poly_idx = bisect_right(self.breaks[0], x) - 1
		poly_shift = self.breaks[0][poly_idx]
		return self.horner(self.coefs[poly_idx], x - poly_shift)

class PPoly():
	
	def mkpp(self, breaks, coefs):	
		d = 1
		sum= len(coefs[0]) * len(coefs)
		dlk = sum
		l = len(breaks)-1
		dl = d*l

		k = dlk/dl+100*(math.pow(2,-52))
		if k<0:
			k = math.ceil(k)
		else:
			k = math.floor(k)

		if k<=0 or (dl*k!=dlk):
			print ("ERROR: MISMATCH PP AND COEF")
			return None
		
		pp = PiecePoly()
		pp.form = 'pp'
		pp.breaks = np.reshape(breaks,(1,l+1))
		pp.coefs = np.reshape(coefs,(dl,k))
		pp.order = k
		pp.dim = d

		return pp