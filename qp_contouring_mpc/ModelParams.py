import numpy as np
import math

class ModelParams():
	def set_info(self, nx, nu, dt_sim):
		self.nx = nx
		self.nu = nu
		self.dt_sim = dt_sim

	def __init__(self, Cm1, Cm2, Cr0, Cr2, Br, Cr, Dr, Bf, Cf, Df, m, Iz, lf, lr):
		self.Cm1 = Cm1
		self.Cm2 = Cm2
		self.Cr0 = Cr0
		self.Cr2 = Cr2
		self.Br = Br
		self.Cr = Cr
		self.Dr = Dr
		self.Bf = Bf
		self.Cf = Cf
		self.Df = Df
		self.m = m
		self.Iz = Iz
		self.lf = lf
		self.lr = lr
		
		self.Tx = np.diag(1 / np.array([3,3,2*math.pi,4,2,7,30]))
		self.invTx = np.linalg.inv(self.Tx)
		self.Tu = np.diag(1 / np.array([1,0.35,6]))
		self.invTu = np.linalg.inv(self.Tu)
