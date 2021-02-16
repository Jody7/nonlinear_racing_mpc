from ModelParams import ModelParams
from SimulateSys import SimulateSys
from DiscreteLinearModel import DiscreteLinearModel

import numpy as np
import math

def start_mpc():
	sim = SimulateSys()
	dslm = DiscreteLinearModel()

	dt_sim = 0.01
	x = np.zeros(7)
	u = np.zeros(3)
	model_params = ModelParams(0.287, 0.0545, 0.0518, 0.0035, 3.38, 1.26, 0.173, 2.57, 1.2, 0.192, 0.041, 27.8e-6, 0.029, 0.033)

	# global control loop
	stop_time = 1.0
	for global_sim_time in np.linspace(0, stop_time, int(stop_time/dt_sim)):
		u[0] = 0.3
		u[1] = -0.1

		dslm.get_discrete_model(x, u, model_params, dt_sim)

		# perform global sim step
		x_solved_state = sim.step(dt_sim, x, u, model_params)
		x = x_solved_state

np.set_printoptions(precision=4, suppress=True)
start_mpc()