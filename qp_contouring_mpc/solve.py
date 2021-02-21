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

def augmentState(sim, x, u, x0, model_params, N):
	x[:,0] = x0
	u[:,0] = u[:,1]

	for i in range(N-1):
		if i != 0:
			x[:,i] = x[:,i+1]
			u[:,i] = u[:,i+1]
	x[:,N-1] = x[:,N]
	u[:,N-1] = u[:,N-1]	# :3 helps readibility

	x[:,N] = sim.step(model_params.dt_sim, x[:,N], u[:,N-1], model_params)

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


def start_mpc():
	global x, u, x0, u0, car_pos_list_x, car_pos_list_y
	sim = SimulateSys()

	dt_sim = 0.01
	x0 = np.zeros(7)
	u0 = np.zeros(3)
	model_params = ModelParams(0.287, 0.0545*0.8, 0.0518*0.8, 0.0035*0.8, 3.38, 1.26, 0.173, 2.57, 1.2, 0.192, 0.041, 27.8e-6, 0.029, 0.033)
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
      0,0,0,theta])
	print(x0)
    # ------------

	mpc_solve = MPC_solve(dlsm)
	N = 10
	x = np.zeros((7,N+1)) # needs to be initialized as copys of x0 eventually
	u = np.zeros((3,N))

	fig, (ax1) = plt.subplots(1,1)
	lines = [(ax1.plot([], [], lw=2)[0])]
	car_pos_list_x = []
	car_pos_list_y = []
	axis_lims_inner = 4
	ax1.set_xlim(-axis_lims_inner, axis_lims_inner)
	ax1.set_ylim(-axis_lims_inner, axis_lims_inner)
	ax1.set_aspect('equal', 'box')

	def solve(solve_idx):
		global x, u, x0, u0
		x, u = augmentState(sim, x, u, x0, model_params, N)
		x, u = mpc_solve.solve_routine(TrackMPC, N, model_params, x, u, x0, u0)

		print("---")
		for i in range(1):
			print("x", x[:,i])
			if i < N:
				print("u", u[:,i])

		# apply linearized mp estimate and solution to ODE to get nonlinear sol
		x0 = sim.step(dt_sim, x[:,0], u[:,0], model_params)
		u0 = u[:,0]
		theta = find_theta([x0[0], x0[1]],track['center'],traj['ppx']['breaks'],trackWidth,0)
		print("theta:", theta)
		x0[model_params.nx-1] = theta

		car_pos_list_x.append(x[0][0])
		car_pos_list_y.append(x[1][0])
		lines[0].set_data(car_pos_list_x, car_pos_list_y)

		return lines


	anim = animation.FuncAnimation(fig, solve, 100000, interval=1, blit=True, repeat=False)
	plt.show()

	# global control loop
	#stop_time = 1.0
	#for global_sim_time in np.linspace(0, stop_time, int(stop_time/dt_sim)):
	#	solve()


np.set_printoptions(precision=4, suppress=True)
load_matlab()
start_mpc()