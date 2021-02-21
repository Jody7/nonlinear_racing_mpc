import math
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
		poly_idx = bisect_right(self.breaks[0], x)
		return self.horner(self.coefs[poly_idx], x)

class PPoly():
	def mkpp(self, breaks,coefs,*args):
		if len(args)==1:
			d = np.transpose(args[0])
		else:
			d = 1
		sum=0
		try:
			for i in range(len(coefs)):
				for j in range(len(coefs[i])):
					sum = sum+len(coefs[i][j])
		except:
			try:
				for i in range(len(coefs)):
					sum = sum+len(coefs[i])
			except:
				sum = len(coefs)

		dlk = sum
		l = len(breaks)-1

		try:
			if len(d) > 1:
				prod = 0
				for i in range(len(d)):
					prod = prod*d[i]
				dl = prod*l
			else:
				dl = d*l
		except:
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