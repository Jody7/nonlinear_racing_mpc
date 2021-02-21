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
