from DiscreteLinearModel import DiscreteLinearModel

import time, math, ctypes
import numpy as np
from PPoly import PPoly
from scipy.linalg import block_diag

ppoly = PPoly()

bounds = np.matrix([
	[-1, 1],
	[-1, 1],
	[-3, 3],
	[0, 0.2],
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

def stage_flatten(stage_idx, flat_mat, dest):
	for elem_idx in range(len(flat_mat)):
		dest[elem_idx + len(flat_mat) * stage_idx] = flat_mat[elem_idx]
def stage_select(orig_mat, sel_opts):
	res = []
	for sel in sel_opts:
		res.append(orig_mat[sel[0]-1, sel[1]-1])
	return res

class MPC_solve():

	def __init__(self, dlsm):
		self.dlsm = dlsm
		self.so_file = "testsolver.exe"
		self.so_load = ctypes.CDLL(self.so_file)
		self.so_load.perform_solve.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_double))

	def solve_routine(self, TrackMPC, N, model_params, Xhor, Uhor, x0, u0):
		t0 = time.time()

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
				stage[i]["Ck"], stage[i]["ug"], stage[i]["lg"] = getInequalityConstraints(borders[:,max(i-1, 0)], model_params)
			stage[i]["Qk"] = generateH(TrackMPC, model_params, Xk, i, N)
			stage[i]["fk"] = generateF(TrackMPC, model_params, Xk, i, N)
			stage[i]["Rk"] = 2 * np.diag([0.01, 0.1, 1e-3])

		T = N - 1

		_Qk = [0] * 100 * (T+2)
		_Rk = [0] * 3 * (T+1)
		_fk = [0] * 4 * (T+2)
		_Ak = [0] * 40 * (T+1)
		_Bk = [0] * 16 * (T+1)
		_gk = [0] * 6 * (T+1)
		_Ck = [0] * 2 * (T+2)
		_ug = [0] * 1 * (T+2)
		_lg = [0] * 1 * (T+2)
		_x = np.asarray((np.diag(np.concatenate((model_params.Tx.diagonal(), model_params.Tu.diagonal()))) @ np.matrix(np.concatenate(([stage[0]["x0"], stage[0]["u0"]]))).T).flatten())[0]

		for stage_idx in range(N):
			stage_flatten(stage_idx, stage[stage_idx]["Rk"].diagonal(), _Rk)
			stage_flatten(stage_idx, stage_select(stage[stage_idx]["Ak"], [(1,1), (1,3), (1,4), (1,5), (1,6), (1,8), (1,9), (2,2), (2,3), (2,4), (2,5), (2,6), (2,8), (2,9), (3,3), (3,4), (3,5), (3,6), (3,8), (3,9), (4,4), (4,5), (4,6), (4,8), (4,9), (5,4), (5,5), (5,6), (5,8), (5,9), (6,4), (6,5), (6,6), (6,8), (6,9), (7,7), (7,10), (8,8), (9,9), (10,10)]), _Ak)
			stage_flatten(stage_idx, stage_select(stage[stage_idx]["Bk"], [(1,1), (1,2), (2,1), (2,2), (3,1), (3,2), (4,1), (4,2), (5,1), (5,2), (6,1), (6,2), (7,3), (8,1), (9,2), (10,3)]), _Bk)
			stage_flatten(stage_idx, stage_select(stage[stage_idx]["gk"], [(1,1), (2,1), (3,1), (4,1), (5,1), (6,1)]), _gk)

		for stage_idx in range(N+1):
			stage_flatten(stage_idx, stage[stage_idx]["Qk"].flatten('C'), _Qk)
			stage_flatten(stage_idx, stage_select(stage[stage_idx]["fk"], [(1,1), (2,1), (7,1), (10,1)]), _fk)
			stage_flatten(stage_idx, [stage[stage_idx]["Ck"][0], stage[stage_idx]["Ck"][1]], _Ck)
			stage_flatten(stage_idx, [stage[stage_idx]["ug"][0, 0]], _ug)
			stage_flatten(stage_idx, [stage[stage_idx]["lg"][0, 0]], _lg)

		t_0 = time.time()
		res = self.so_load.perform_solve(T,
			(ctypes.c_double*len(_Qk))(*_Qk),
			(ctypes.c_double*len(_Rk))(*_Rk),
			(ctypes.c_double*len(_fk))(*_fk),
			(ctypes.c_double*len(_Ak))(*_Ak),
			(ctypes.c_double*len(_Bk))(*_Bk),
			(ctypes.c_double*len(_gk))(*_gk),
			(ctypes.c_double*len(_Ck))(*_Ck),
			(ctypes.c_double*len(_ug))(*_ug),
			(ctypes.c_double*len(_lg))(*_lg),
			(ctypes.c_double*len(_x))(*_x)
		)
		#print("solver hz = " , 1/(time.time() - t_0))

		x_solved_list = []
		for i in range(N+1):
			x_solved_list.append([res[i][0], res[i][1], res[i][2], res[i][3], res[i][4], res[i][5], res[i][6], res[i][7], res[i][8], res[i][9]])
		x_solved_mat = np.array(x_solved_list).transpose()

		x_value = model_params.invTx @ x_solved_mat[0:model_params.nx,]
		u_value = model_params.invTu @ x_solved_mat[model_params.nx:model_params.nx+model_params.nu,1:]

		#matprint(x_value, ".2f")

		return x_value, u_value, borders