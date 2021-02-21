from DiscreteLinearModel import DiscreteLinearModel

import cvxpy, time
import numpy as np
from PPoly import PPoly

ppoly = PPoly()

def getEqualityConstraints(Xk, Uk, model_params, dlsm):
	nx = model_params.nx
	nu = model_params.nu

	Ad, Bd, gd = dlsm.get_discrete_model(Xk, Uk, model_params, model_params.dt_sim)
	# no input/state scaling
	return Ad, Bd, gd

def getInequalityConstraints(border, model_params):
	nx = model_params.nx
	nu = model_params.nu

	x1 = border[0][0]
	y1 = border[1][0]
	x2 = border[2][0]
	y2 = border[3][0]

	numer = -(x2-x1)
	denom = (y2-y1)
	dbmax = max(numer*x1-denom*y1,numer*x2-denom*y2)
	dbmin = min(numer*x1-denom*y1,numer*x2-denom*y2)

	Ck = np.zeros(nx+nu)
	Ck[0] = numer
	Ck[1] = -denom
	ug = dbmax
	lg = dbmin

	return Ck, ug, lg


class MPC_solve():

	def __init__(self, dlsm):
		print(cvxpy.installed_solvers())
		self.dlsm = dlsm

	def solve_routine(self, TrackMPC, N, model_params, Xhor, Uhor, x0, u0):
		theta_phys = []
		for i in range(N):
			theta_virt = Xhor[model_params.nx - 1][i] % TrackMPC['tl']
			theta_phys.append(theta_virt)

		TrackLeftx = []
		TrackLefty = []
		TrackRightx = []
		TrackRighty = []

		for theta_phy in theta_phys:
			TrackLeftx.append(ppoly.mkpp(TrackMPC['borders']['pplx']['breaks'], TrackMPC['borders']['pplx']['coefs']).eval(theta_phy))
			TrackLefty.append(ppoly.mkpp(TrackMPC['borders']['pply']['breaks'], TrackMPC['borders']['pply']['coefs']).eval(theta_phy))
			TrackRightx.append(ppoly.mkpp(TrackMPC['borders']['pprx']['breaks'], TrackMPC['borders']['pprx']['coefs']).eval(theta_phy))
			TrackRighty.append(ppoly.mkpp(TrackMPC['borders']['ppry']['breaks'], TrackMPC['borders']['ppry']['coefs']).eval(theta_phy))

		borders = np.matrix([TrackLeftx, TrackLefty, TrackRightx, TrackRighty])

		x = cvxpy.Variable((model_params.nx, N+1))
		u = cvxpy.Variable((model_params.nu, N))

		stage = [{}] * N
		stage[0]["x0"] = x0
		stage[0]["u0"] = u0

		for i in range(0, N):
			Xk = Xhor[:,i];
			Uk = Uhor[:,i];
			stage[0]["Ak"], stage[0]["Bk"], stage[0]["gk"] = getEqualityConstraints(Xk, Uk, model_params, self.dlsm)
			stage[0]["Ck"], stage[0]["ug"], stage[0]["lg"] = getInequalityConstraints(borders[:,max(i, 1)], model_params)

		cost = 0.0
		constr = []
		constr += [x[:,0] == stage[0]["x0"]]
		constr += [u[:,0] == stage[0]["u0"]]
		for stage_idx in range(len(stage)):
			constr += [x[:,stage_idx+1] == stage[i]["Ak"]@x[:,stage_idx] + stage[i]["Bk"]@u[:,stage_idx]]

			constr_u_idx = min(stage_idx+1, len(stage)-1)
			constr += [u[0][constr_u_idx] - u[0][stage_idx] <= 3 * model_params.dt_sim]
			constr += [u[0][constr_u_idx] - u[0][stage_idx] >= -3 * model_params.dt_sim]
			constr += [u[1][constr_u_idx] - u[1][stage_idx] <= 0.1 * model_params.dt_sim]
			constr += [u[1][constr_u_idx] - u[1][stage_idx] >= -0.1 * model_params.dt_sim]

			constr += [u[0][stage_idx] >= -1.0]
			constr += [u[0][stage_idx] <= 1.0]
			constr += [u[1][stage_idx] >= -0.1]
			constr += [u[1][stage_idx] <= 0.1]

			#cost += (1 - x[0][stage_idx])**2
			#cost += (1 - x[1][stage_idx])**2
			cost += (1.0 - u[0][stage_idx])**2
			cost += (-0.2 - u[1][stage_idx])**2

		problem = cvxpy.Problem(cvxpy.Minimize(cost), constr)
		starttime = time.time()
		problem.solve(solver=cvxpy.GUROBI, warm_start=True, verbose=False)
		print("python resolved result: ", time.time()-starttime)

		return x.value, u.value