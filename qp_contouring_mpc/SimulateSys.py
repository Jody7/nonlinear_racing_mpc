import numpy as np
import math
from scipy.integrate import solve_ivp

class SimulateSys():
	def calc_dynamics(self, t, x, u, ModelParams):
		Cm1=ModelParams.Cm1
		Cm2=ModelParams.Cm2
		Cr0=ModelParams.Cr0
		Cr2=ModelParams.Cr2
				
		B_r = ModelParams.Br
		C_r = ModelParams.Cr
		D_r = ModelParams.Dr
		B_f = ModelParams.Bf
		C_f = ModelParams.Cf
		D_f = ModelParams.Df
		
		m = ModelParams.m
		Iz = ModelParams.Iz
		l_f = ModelParams.lf
		l_r = ModelParams.lr

		phi   =x[2]
		v_x   =x[3]
		v_y   =x[4]
		omega =x[5]
		D     =u[0]
		delta =u[1]
		
		alpha_f = -math.atan2(l_f*omega + v_y, abs(v_x)) + delta
		alpha_r =  math.atan2(l_r*omega - v_y, abs(v_x))

		F_fy = D_f*math.sin(C_f*math.atan(B_f*alpha_f))
		F_ry = D_r*math.sin(C_r*math.atan(B_r*alpha_r))

		F_rx = (Cm1*D-Cm2*D*v_x-Cr0-Cr2*v_x**2)

		xdot=np.array([v_x*math.cos(phi) - v_y*math.sin(phi),
			v_y*math.cos(phi) + v_x*math.sin(phi),
			omega,
			1/m*(F_rx - F_fy*math.sin(delta) + m*v_y*omega),
			1/m*(F_ry + F_fy*math.cos(delta) - m*v_x*omega),
			1/Iz*(F_fy*l_f*math.cos(delta)- F_ry*l_r),
			u[2]])

		return xdot

	def step(self, dt, x, u, ModelParams):
		sol = solve_ivp(self.calc_dynamics, (0, dt), x, args=(u, ModelParams), method = 'RK45')
		return sol.y[:, -1]