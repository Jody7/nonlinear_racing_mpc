from ModelParams import ModelParams
from SimulateSys import SimulateSys
from MPC_solve import  MPC_solve
from DiscreteLinearModel import DiscreteLinearModel
from MatLoader import MatLoader
from PPoly import PPoly

import numpy as np
import math, time

import matplotlib.pyplot as plt
from matplotlib import animation

matloader = MatLoader()
ppoly = PPoly()

def augmentState(sim, x, u, x0, model_params, N, tl):
	x[:,0] = x0
	u[:,0] = u[:,1]

	for i in range(N-1):
		if i != 0:
			x[:,i] = x[:,i+1]
			u[:,i] = u[:,i+1]
	x[:,N-1] = x[:,N]
	u[:,N-1] = u[:,N-1]	# :3 helps readibility

	x[:,N] = sim.step(model_params.dt_sim, x[:,N], u[:,N-1], model_params)

	if x[2,0] - x[2,1] > math.pi:
		x[2,1:] = x[2,1:] + 2*math.pi
	if x[2,0] - x[2,1] < -math.pi:
		x[2,1:] = x[2,1:] - 2*math.pi
	if x[-1,0] - x[-1,1] < -0.75*tl:
		x[-1,1:] = x[-1,1:] - tl

	return x, u

def load_matlab():
	track_data = matloader.loadmat('track.mat')['track']
	traj_data = matloader.loadmat('traj.mat')['traj']
	border_data = matloader.loadmat('borders.mat')['borders']
	track_len = traj_data['ppy']['breaks'][-1]

	return {'traj' : traj_data, 'borders' : border_data, 'tl' : track_len, 'track' : track_data}

def find_theta(currentPosition,trackCenter,traj_breaks,trackWidth,last_closestIdx):
	posX = currentPosition[0]
	posY = currentPosition[1]
	smallest_dist = 9999
	smallest_idx = 0
	for i in range(len(trackCenter[0])):
		search_idx = (i + last_closestIdx - 20) % len(trackCenter[0])
		search_dist = (trackCenter[0][search_idx] - posX)**2 + (trackCenter[1][search_idx] - posY)**2
		# TODO: optimization with locallity search stuff, can either detect minima or do simple dist bounds
		if search_dist < smallest_dist:
			smallest_dist = search_dist
			smallest_idx = search_idx
	previous_idx = (smallest_idx - 1) % len(trackCenter[0])
	cosinus = np.dot(np.array([posX, posY]) - trackCenter[:,smallest_idx], trackCenter[:,previous_idx] - trackCenter[:,smallest_idx])
	minIndex2 = 0
	minIndex = smallest_idx

	if (cosinus > 0):
		minIndex2 = previous_idx
	else:
		minIndex2 = smallest_idx
		minIndex = (smallest_idx+1) % len(trackCenter[0])
	
	if abs(smallest_dist) != 0:
		cosinus = np.dot(np.array([posX, posY]) - trackCenter[:,minIndex2], trackCenter[:,minIndex] - trackCenter[:,minIndex2]) / (np.linalg.norm(np.array([posX, posY]) - trackCenter[:,minIndex2])*np.linalg.norm(trackCenter[:,minIndex] - trackCenter[:,minIndex2]))
	else:
		cosinus = 0
	theta = traj_breaks[minIndex2]
	theta = theta + cosinus*np.linalg.norm(np.array([posX, posY]) - trackCenter[:,minIndex2], 2)
	return theta 

def unwrap_x0(x0):
	if x0[2] > math.pi:
		x0[2] = x0[2] - 2*math.pi
	if x0[2] < -math.pi:
		x0[2] = x0[2] + 2*math.pi
	return x0

def start_mpc():
	global x, u, x0, u0, car_pos_list_x, car_pos_list_y
	sim = SimulateSys()

	dt_sim = 0.025
	x0 = np.zeros(7)
	u0 = np.zeros(3)
								#  Cm1, Cm2,   Cr0,     Cr2,     Br,          Cr , Dr,    Bf,    Cf,   Df,    m,    Iz,       lf,     lr
	model_params = ModelParams(0.287*1.5, 0.0545*3, 0.0518, 0.00035, 3.3852*0.5,  0.8, 0.1737, 2.579, 1.2, 0.192, 0.041, 27.8e-6, 0.029, 0.033)
	model_params.set_info(7, 3, dt_sim)
	dlsm = DiscreteLinearModel()

	# ----------- SIM LOAD -----------

	TrackMPC = load_matlab()
	track = TrackMPC['track']
	traj = TrackMPC['traj']
	trackWidth = 0.28
	startIdx = 0

	theta = find_theta([track['center'][0][startIdx],track['center'][1][startIdx]],track['center'],traj['ppx']['breaks'],trackWidth,startIdx)

	x0 = np.array([track['center'][0][startIdx],track['center'][1][startIdx],
      	math.atan2(ppoly.mkpp(traj['dppy']['breaks'], traj['dppy']['coefs']).eval(theta),
		ppoly.mkpp(traj['dppx']['breaks'], traj['dppx']['coefs']).eval(theta)),
      0.0,0,0,theta])
    # ------------

	mpc_solve = MPC_solve(dlsm)
	N = 25
	x = np.zeros((7,N+1))
	for i in range(N+1):
		x[:,i] = x0
	u = np.zeros((3,N))

	fig, (ax1) = plt.subplots(1,1)
	lines = [(ax1.plot([], [], lw=3)[0]), (ax1.plot([], [], linestyle='dashed')[0]), (ax1.plot([], [], linestyle='dashed')[0]), (ax1.plot([], [], linestyle='dotted')[0]), (ax1.plot([], [], linestyle='dotted')[0]), (ax1.plot([], [], lw=2, color='pink')[0]) , (ax1.plot([], [], lw=2, color='pink')[0])]
	lines_len = len(lines)
	lines.append(None)

	car_pos_list_x = []
	car_pos_list_y = []
	axis_lims_inner = 2
	ax1.set_xlim(-axis_lims_inner, axis_lims_inner)
	ax1.set_ylim(-axis_lims_inner, axis_lims_inner)
	ax1.set_aspect('equal', 'box')

	lines[1].set_data(TrackMPC['track']['inner'][0], TrackMPC['track']['inner'][1])
	lines[2].set_data(TrackMPC['track']['outer'][0], TrackMPC['track']['outer'][1])

	def solve(solve_idx):
		global x, u, x0, u0
		x, u = augmentState(sim, x, u, x0, model_params, N, TrackMPC['tl'])
		borders = None

		sqp_damping_x = 0.7
		sqp_damping_u = 0.7
		for sqp_idx in range(2):
			x_solve, u_solve, borders = mpc_solve.solve_routine(TrackMPC, N, model_params, x, u, x0, u0)
			x = sqp_damping_x*x + (1-sqp_damping_x)*x_solve
			u = sqp_damping_u*u + (1-sqp_damping_u)*u_solve

		# apply linearized mp estimate and solution to ODE to get nonlinear sol
		x0 = sim.step(dt_sim, x[:,0], u[:,0], model_params)
		x0 = unwrap_x0(x0)
		u0 = u[:,0]
		theta = find_theta([x0[0], x0[1]],track['center'],traj['ppx']['breaks'],trackWidth,0)
		#print("theta:", theta, "u0 in solve:" ,u0, x0[2])
		x0[model_params.nx-1] = theta

		if len(car_pos_list_x) > 100:
			for i in range(10):
				car_pos_list_x.pop(0)
				car_pos_list_y.pop(0)

		car_pos_list_x.append(x[0][0])
		car_pos_list_y.append(x[1][0])
		lines[0].set_data(car_pos_list_x, car_pos_list_y)
		lines[3].set_data(x[0], x[1])
		
		plot_traj_data_x = []
		plot_traj_data_y = []
		for i in range(len(x[6])):
			plot_traj_data_x.append(ppoly.mkpp(TrackMPC['traj']['ppx']['breaks'], TrackMPC['traj']['ppx']['coefs']).eval(x[6][i]))
			plot_traj_data_y.append(ppoly.mkpp(TrackMPC['traj']['ppy']['breaks'], TrackMPC['traj']['ppy']['coefs']).eval(x[6][i]))
		lines[4].set_data(plot_traj_data_x, plot_traj_data_y)
		lines[5].set_data(borders[0], borders[1]) # Left Side
		lines[6].set_data(borders[2], borders[3]) # Right Side

		car_arrow_size = 0.1
		car_arrow = plt.Arrow(x[0][0],x[1][0],car_arrow_size*math.cos(x[2][0]),car_arrow_size*math.sin(x[2][0]), color = 'r', width=0.3)
		ax1.add_patch(car_arrow)
		lines[lines_len] = car_arrow

		return lines

	anim = animation.FuncAnimation(fig, solve, 100000, interval=1, blit=True, repeat=False)
	plt.show()


np.set_printoptions(precision=4, suppress=True)
np.show_config()
load_matlab()
start_mpc()