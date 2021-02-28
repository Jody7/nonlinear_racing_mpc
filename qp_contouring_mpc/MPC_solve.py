from DiscreteLinearModel import DiscreteLinearModel

import cvxpy, time, math
import numpy as np
from PPoly import PPoly
from scipy.linalg import block_diag

ppoly = PPoly()

bounds = np.matrix([
	[-1, 1],
	[-1, 1],
	[-500, 500],
	[0, 1],
	[-1, 1],
	[-1, 1],
	[0, 1],
	[-0.1, 1],
	[-1, 1],
	[0, 1],
	[-1, 1],
	[-1, 1],
	[-5, 5],
	])

def matprint(mat, fmt="g"):
	col_maxes = [max([len(("{:"+fmt+"}").format(x)) for x in col]) for col in mat.T]
	for x in mat:
		for i, y in enumerate(x):
			print(("{:"+str(col_maxes[i])+fmt+"}").format(y), end="  ")
		print("")

def getdphivirt_dtheta(theta_virt, TrackMPC):
	dxdth = ppoly.mkpp(TrackMPC['traj']['dppx']['breaks'], TrackMPC['traj']['dppx']['coefs']).eval(theta_virt)
	dydth = ppoly.mkpp(TrackMPC['traj']['dppy']['breaks'], TrackMPC['traj']['dppy']['coefs']).eval(theta_virt)
	d2xdth2 = ppoly.mkpp(TrackMPC['traj']['ddppx']['breaks'], TrackMPC['traj']['ddppx']['coefs']).eval(theta_virt)
	d2ydth2 = ppoly.mkpp(TrackMPC['traj']['ddppy']['breaks'], TrackMPC['traj']['ddppy']['coefs']).eval(theta_virt)

	numer=dxdth*d2ydth2 - dydth*d2xdth2
	denom=dxdth**2 + dydth**2

	return numer / denom

def getderror_dtheta(TrackMPC, theta_virt, x_phys, y_phys):
	dxvirt_dtheta = ppoly.mkpp(TrackMPC['traj']['dppx']['breaks'], TrackMPC['traj']['dppx']['coefs']).eval(theta_virt)
	dyvirt_dtheta = ppoly.mkpp(TrackMPC['traj']['dppy']['breaks'], TrackMPC['traj']['dppy']['coefs']).eval(theta_virt)

	phi_virt = math.atan2(dyvirt_dtheta, dxvirt_dtheta)

	x_virt = ppoly.mkpp(TrackMPC['traj']['ppx']['breaks'], TrackMPC['traj']['ppx']['coefs']).eval(theta_virt)
	y_virt = ppoly.mkpp(TrackMPC['traj']['ppy']['breaks'], TrackMPC['traj']['ppy']['coefs']).eval(theta_virt)

	Dx = x_phys - x_virt
	Dy = y_phys - y_virt

	dphivirt_dtheta = getdphivirt_dtheta(theta_virt, TrackMPC)
	cos_phi_virt = math.cos(phi_virt)
	sin_phi_virt = math.sin(phi_virt)

	tmp1 = np.array([dphivirt_dtheta, 1])
	tmp2 = np.array([cos_phi_virt, sin_phi_virt])

	MC = np.matrix(
		[[Dx, Dy],
		[dyvirt_dtheta, -dxvirt_dtheta]])
	ML = np.matrix(
		[[-Dy, Dx],
		[dxvirt_dtheta, dyvirt_dtheta]])

	deC_dtheta = tmp1 @ MC @ tmp2
	deL_dtheta = tmp1 @ ML @ tmp2

	return deC_dtheta.item(), deL_dtheta.item(), cos_phi_virt, sin_phi_virt


def getErrorGradient(TrackMPC, theta_virt, model_params, x_phys, y_phys):
	deC_dtheta, deL_dtheta, cos_phi_virt, sin_phi_virt = getderror_dtheta(TrackMPC, theta_virt, x_phys, y_phys)

	grad_eC = np.array([sin_phi_virt, -cos_phi_virt, 0, 0, 0, 0, deC_dtheta])
	grad_eL = np.array([-cos_phi_virt, -sin_phi_virt, 0, 0, 0, 0, deL_dtheta])

	return grad_eC, grad_eL


def generateQtilde(TrackMPC, model_params, Xk, i, N):
	if i == N:
		Q = np.diag([10 * 0.1, 1000])
	else:
		Q = np.diag([0.1, 1000])

	theta_virt = Xk[-1] % TrackMPC['traj']['ppx']['breaks'][-1]
	grad_eC, grad_eL = getErrorGradient(TrackMPC, theta_virt, model_params, Xk[0], Xk[1])
	errorgrad = np.matrix(
		[grad_eC,
		grad_eL]) # writing this as vstack probably is faster
	Qtilde = errorgrad.transpose() @ Q @ errorgrad

	return Qtilde


def generateH(TrackMPC, model_params, Xk, i, N):
	Qtilde = generateQtilde(TrackMPC, model_params, Xk, i, N)
	if i == N:
		Qtilde[5, 5] = 10 * 1e-05
	else:
		Qtilde[5, 5] = 1e-05

	Qtilde = 0.5 * (Qtilde + Qtilde.transpose())
	Qk = 2*block_diag(Qtilde, np.diag([1e-4, 1e-4, 1e-4]))
	Qk = block_diag(model_params.invTx, model_params.invTu) @ Qk @ block_diag(model_params.invTx, model_params.invTu) + 1e-12*np.eye(10)
	return Qk

def getErrors(TrackMPC, theta_virt, x_phys, y_phys):
	dxdth = ppoly.mkpp(TrackMPC['traj']['dppx']['breaks'], TrackMPC['traj']['dppx']['coefs']).eval(theta_virt)
	dydth = ppoly.mkpp(TrackMPC['traj']['dppy']['breaks'], TrackMPC['traj']['dppy']['coefs']).eval(theta_virt)

	x_virt = ppoly.mkpp(TrackMPC['traj']['ppx']['breaks'], TrackMPC['traj']['ppx']['coefs']).eval(theta_virt)
	y_virt = ppoly.mkpp(TrackMPC['traj']['ppy']['breaks'], TrackMPC['traj']['ppy']['coefs']).eval(theta_virt)

	print(theta_virt,dxdth,dydth,x_virt,y_virt);

	phi_virt = math.atan2(dydth, dxdth)
	sin_phi_virt = math.sin(phi_virt)
	cos_phi_virt = math.cos(phi_virt)

	eC = -sin_phi_virt*(x_virt - x_phys) + cos_phi_virt*(y_virt - y_phys)
	eL =  cos_phi_virt*(x_virt - x_phys) + sin_phi_virt*(y_virt - y_phys)

	return eC, eL

def generateF(TrackMPC, model_params, Xk, i, N):
	x_phys = Xk[0]
	y_phys = Xk[1]

	theta_virt = Xk[-1] % TrackMPC['traj']['ppx']['breaks'][-1]
	eC, eL = getErrors(TrackMPC, theta_virt, x_phys, y_phys)
	e = np.vstack((eC, eL))
	grad_eC, grad_eL = getErrorGradient(TrackMPC, theta_virt, model_params, x_phys, y_phys)
	grad_e = np.matrix(
		[grad_eC,
		grad_eL])

	print(e)

	if i == N:
		Q = np.diag([10 * 0.1, 1000])
	else:
		Q = np.diag([0.1, 1000])

	fx = (2 * e.transpose() @ Q @ grad_e) - (2 * Xk.transpose() @ grad_e.transpose() @ Q @ grad_e)
	fT = np.hstack((fx, [np.zeros(model_params.nu - 1)], [[-0.02]]))
	f = block_diag(model_params.invTx, model_params.invTu) @ fT.transpose()

	return f
		
def getEqualityConstraints(Xk, Uk, model_params, dlsm):
	nx = model_params.nx
	nu = model_params.nu

	Ad, Bd, gd = dlsm.get_discrete_model(Xk, Uk, model_params, model_params.dt_sim)

	Ak = np.hstack((model_params.Tx @ Ad @ model_params.invTx, model_params.Tx @ Bd @ model_params.invTu))
	Ak = np.vstack((Ak, np.hstack((np.zeros((model_params.nu, model_params.nx)), np.eye(3)))))
	Bk = np.vstack((model_params.Tx @ Bd @ model_params.invTu, np.eye(model_params.nu)))
	gk = np.vstack((model_params.Tx @ gd, np.matrix(np.zeros(model_params.nu)).transpose()))

	return Ak, Bk, gk

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
	Ck = Ck @ np.diag(np.concatenate((model_params.invTx.diagonal(), model_params.invTu.diagonal())))

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

		x = cvxpy.Variable((model_params.nx+model_params.nu, N+1))
		u = cvxpy.Variable((model_params.nu, N))

		stage = {}
		for i in range (0, N+1):
			stage[i] = {}

		stage[0]["x0"] = x0
		stage[0]["u0"] = u0

		for i in range(0, N+1):
			Xk = Xhor[:,i];
			if i < N:
				Uk = Uhor[:,i];
				stage[i]["Ak"], stage[i]["Bk"], stage[i]["gk"] = getEqualityConstraints(Xk, Uk, model_params, self.dlsm)
				stage[i]["Ck"], stage[i]["ug"], stage[i]["lg"] = getInequalityConstraints(borders[:,max(i-1, 0)], model_params)
			else:
				stage[i]["Ck"], stage[i]["ug"], stage[i]["lg"] = getInequalityConstraints(borders[:,max(i-2, 0)], model_params)

			stage[i]["Qk"] = generateH(TrackMPC, model_params, Xk, i, N)
			stage[i]["fk"] = generateF(TrackMPC, model_params, Xk, i, N)
			stage[i]["Rk"] = 2 * np.diag([0.01, 1, 1e-3])

		cost = 0.0
		constr = []
		constr += [x[:,0] == cvxpy.reshape(np.diag(np.concatenate((model_params.Tx.diagonal(), model_params.Tu.diagonal()))) @ np.matrix(np.concatenate(([stage[0]["x0"], stage[0]["u0"]]))).T, (10, ))]

		for stage_idx in range(N):
			#matprint(stage[stage_idx]["Qk"], ".6f")
			#matprint(stage[stage_idx]["Rk"], ".6f")
			print("fk", stage_idx, stage[stage_idx]["fk"][0])
			#matprint(stage[stage_idx]["Ak"], ".6f")
			#matprint(stage[stage_idx]["Bk"], ".6f")
			#print(stage[stage_idx]["gk"])
			#print(stage[stage_idx]["Ck"])
			#print(stage[stage_idx]["ug"])
			#print(stage[stage_idx]["lg"])
			#print("---", stage_idx)

			constr += [cvxpy.reshape(x[:,stage_idx+1], (10, 1)) == (cvxpy.reshape(stage[stage_idx]["Ak"] @ x[:,stage_idx], (10, 1)) + cvxpy.reshape(stage[stage_idx]["Bk"] @ u[:,stage_idx], (10, 1)) + stage[stage_idx]["gk"])]
			constr += [stage[stage_idx]["Ck"] @ x[:,stage_idx] <= stage[stage_idx]["ug"]]
			constr += [stage[stage_idx]["Ck"] @ x[:,stage_idx] >= stage[stage_idx]["lg"]]

			constr += [cvxpy.reshape(x[:,stage_idx], (model_params.nx+model_params.nu, 1)) >= bounds[:model_params.nx+model_params.nu,0]]
			constr += [cvxpy.reshape(u[:,stage_idx], (model_params.nu, 1)) >= bounds[model_params.nx+model_params.nu:model_params.nx+2*model_params.nu,0]]
			constr += [cvxpy.reshape(x[:,stage_idx], (model_params.nx+model_params.nu, 1)) <= bounds[:model_params.nx+model_params.nu,1]]
			constr += [cvxpy.reshape(u[:,stage_idx], (model_params.nu, 1)) <= bounds[model_params.nx+model_params.nu:model_params.nx+2*model_params.nu,1]]

			cost += 0.5*cvxpy.quad_form(x[:,stage_idx], stage[stage_idx]["Qk"])
			cost += cvxpy.quad_form(u[:,stage_idx], stage[stage_idx]["Rk"])
			cost += (stage[stage_idx]["fk"].transpose() @ x[:,stage_idx])


		# --- AFTERWARDS DO TERMINAL STATE CONSTRAINT AND COST
		stage_idx = N
		cost += 0.5*cvxpy.quad_form(x[:,stage_idx], stage[stage_idx]["Qk"])
		cost += (stage[stage_idx]["fk"].transpose() @ x[:,stage_idx])
		constr += [stage[stage_idx]["Ck"] @ x[:,stage_idx] <= stage[stage_idx]["ug"]]
		constr += [stage[stage_idx]["Ck"] @ x[:,stage_idx] >= stage[stage_idx]["lg"]]
		constr += [cvxpy.reshape(x[:,stage_idx], (model_params.nx+model_params.nu, 1)) >= bounds[:model_params.nx+model_params.nu,0]]
		constr += [cvxpy.reshape(x[:,stage_idx], (model_params.nx+model_params.nu, 1)) <= bounds[:model_params.nx+model_params.nu,1]]

		problem = cvxpy.Problem(cvxpy.Minimize(cost), constr)
		problem.solve(solver=cvxpy.GUROBI, warm_start=True, verbose=False)

		x_value = Xhor
		u_value = Uhor

		print(problem.status)
		if problem.status == "optimal":
			x_value = model_params.invTx @ x.value[0:model_params.nx,]
			u_value = model_params.invTu @ x.value[model_params.nx:model_params.nx+model_params.nu,1:]

			matprint(x_value, ".3f")
			print("==")
			matprint(u_value, ".3f")
			#print("==")
			#matprint(u.value, ".2f")
			print("----------")

		return x_value, u_value