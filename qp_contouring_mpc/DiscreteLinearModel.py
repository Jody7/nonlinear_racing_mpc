import numpy as np
from scipy import signal
import math

class DiscreteLinearModel():
	def get_discrete_model(self, Xbar_k, Ubar_k, ModelParams, Ts):
		Cm1=ModelParams.Cm1
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

		phi     = Xbar_k[2]
		v_x     = Xbar_k[3]
		v_y     = Xbar_k[4]
		omega   = Xbar_k[5]
		D       = Ubar_k[0]
		delta   = Ubar_k[1]
		v_theta = Ubar_k[2]

		if (v_x < 0.5):
			v_x = v_x
			v_y = 0
			omega = 0
			delta = 0
			if(v_x < 0.3):
				v_x = 0.3

		alpha_f = -math.atan2(l_f*omega + v_y,v_x)+delta
		alpha_r =  math.atan2(l_r*omega - v_y,v_x)

		F_fy = D_f*math.sin(C_f*math.atan(B_f*alpha_f))
		F_ry = D_r*math.sin(C_r*math.atan(B_r*alpha_r))

		F_rx = (Cm1*D-Cm2*D*v_x-Cr0-Cr2*v_x**2)

		f = np.row_stack([v_x*math.cos(phi) - v_y*math.sin(phi),
				v_y*math.cos(phi) + v_x*math.sin(phi),
				omega,
				1/m*(F_rx - F_fy*math.sin(delta) + m*v_y*omega),
				1/m*(F_ry + F_fy*math.cos(delta) - m*v_x*omega),
				1/Iz*(F_fy*l_f*math.cos(delta)- F_ry*l_r)])

		# ---------- calculate derivs for discretization

		dFrx_dvx = -Cm2*D - 2*Cr2*v_x;
		dFrx_dD  = Cm1 - Cm2*v_x;

		dFry_dvx = ((B_r*C_r*D_r*math.cos(C_r*math.atan(B_r*alpha_r)))/(1+B_r**2*alpha_r**2)) \
				 *(-(l_r*omega - v_y)/((-l_r*omega + v_y)**2+v_x**2));
			 
		dFry_dvy = ((B_r*C_r*D_r*math.cos(C_r*math.atan(B_r*alpha_r)))/(1+B_r**2*alpha_r**2)) \
				 *((-v_x)/((-l_r*omega + v_y)**2+v_x**2));
			 
		dFry_domega = ((B_r*C_r*D_r*math.cos(C_r*math.atan(B_r*alpha_r)))/(1+B_r**2*alpha_r**2)) \
				 *((l_r*v_x)/((-l_r*omega + v_y)**2+v_x**2));

		dFfy_dvx =     (B_f*C_f*D_f*math.cos(C_f*math.atan(B_f*alpha_f)))/(1+B_f**2*alpha_f**2) \
					 *((l_f*omega + v_y)/((l_f*omega + v_y)**2+v_x**2));
				 
		dFfy_dvy =     (B_f*C_f*D_f*math.cos(C_f*math.atan(B_f*alpha_f)))/(1+B_f**2*alpha_f**2) \
					 *(-v_x/((l_f*omega + v_y)**2+v_x**2));
				 
		dFfy_domega =    (B_f*C_f*D_f*math.cos(C_f*math.atan(B_f*alpha_f)))/(1+B_f**2*alpha_f**2) \
					 *((-l_f*v_x)/((l_f*omega + v_y)**2+v_x**2));
				 
		dFfy_ddelta =  (B_f*C_f*D_f*math.cos(C_f*math.atan(B_f*alpha_f)))/(1+B_f**2*alpha_f**2); 

		df1_dphi = -v_x*math.sin(phi) - v_y*math.cos(phi);
		df1_dvx  = math.cos(phi);
		df1_dvy  = -math.sin(phi);

		df2_dphi = -v_y*math.sin(phi) + v_x*math.cos(phi);
		df2_dvx  = math.sin(phi);
		df2_dvy  = math.cos(phi);

		df3_domega = 1;

		df4_dvx     = 1/m*(dFrx_dvx - dFfy_dvx*math.sin(delta));
		df4_dvy     = 1/m*(           - dFfy_dvy*math.sin(delta)     + m*omega);
		df4_domega = 1/m*(           - dFfy_domega*math.sin(delta) + m*v_y);
		df4_dD      = 1/m*     dFrx_dD;
		df4_ddelta  = 1/m*(           - dFfy_ddelta*math.sin(delta)  - F_fy*math.cos(delta));

		df5_dvx     = 1/m*(dFry_dvx     + dFfy_dvx*math.cos(delta)     - m*omega);
		df5_dvy     = 1/m*(dFry_dvy     + dFfy_dvy*math.cos(delta));   
		df5_domega = 1/m*(dFry_domega + dFfy_domega*math.cos(delta) - m*v_x);
		df5_ddelta  = 1/m*(               dFfy_ddelta*math.cos(delta)  - F_fy*math.sin(delta));

		df6_dvx     = 1/Iz*(dFfy_dvx*l_f*math.cos(delta)     - dFry_dvx*l_r);
		df6_dvy     = 1/Iz*(dFfy_dvy*l_f*math.cos(delta)     - dFry_dvy*l_r);
		df6_domega = 1/Iz*(dFfy_domega*l_f*math.cos(delta) - dFry_domega*l_r);
		df6_ddelta  = 1/Iz*(dFfy_ddelta*l_f*math.cos(delta)  - F_fy*l_f*math.sin(delta));

		# jaaaacoooobi

		Ac = np.row_stack(
			[[0, 0, df1_dphi, df1_dvx, df1_dvy, 0], 
			[0, 0, df2_dphi, df2_dvx, df2_dvy, 0],
			[0, 0, 0, 0, 0, df3_domega],
			[0, 0, 0, df4_dvx, df4_dvy, df4_domega],
			[0, 0, 0, df5_dvx, df5_dvy, df5_domega],
			[0, 0, 0, df6_dvx, df6_dvy, df6_domega]])

		Bc = np.row_stack(
			[[0, 0],
			[0, 0],
			[0, 0],
			[df4_dD, df4_ddelta],
			[0, df5_ddelta],
			[0, df6_ddelta]])

		sx = len(Xbar_k) - 1
		su = len(Ubar_k) - 1

		# need to verify correctness of this
		gc = f - (Ac @ np.row_stack(Xbar_k[0:sx])) - (Bc @ np.row_stack(Ubar_k[0:su]))
		Bc_aug = np.hstack((Bc, gc))

		# super ghetto but should be fine
		ss = signal.StateSpace(Ac, Bc_aug, np.zeros_like(Ac), np.zeros_like(Bc_aug))
		ss = ss.to_discrete(0.1, method='zoh') # zoh should be fine and (i think) this implemention uses matrix exponentiation to do it 

		Ad_ss = np.matrix(ss.A)
		Ad = np.pad(Ad_ss[:,:sx], (0, 1), mode='constant')
		Ad[sx][sx] = 1.0

		Bd_ss = np.matrix(ss.B)
		Bd = np.pad(Bd_ss[:,:su], (0, 1), mode='constant')
		Bd[sx][su] = Ts

		gd = np.pad(Bd_ss[:,-1], ((0, 1), (0, 0)), mode='constant')

		return Ad, Bd, gd