from pyomo.environ import *
from pyomo.dae import *
import pyomo.contrib.fbbt.interval as interval

import matplotlib.pyplot as plt
from matplotlib import animation

import scipy.io

import numpy as np
import math, time

track_mat = scipy.io.loadmat('trackMobil.mat')

solver = SolverFactory('ipopt')
solver.options['tol'] = 1e-4
solver.options['warm_start_init_point'] = 'yes'
solver.options['warm_start_bound_push'] = 1e-2
solver.options['warm_start_mult_bound_push'] = 1e-2
solver.options['expect_infeasible_problem'] = 'no' # don't care if its infeasible
solver.options['mu_init'] = 1e-2

epsilon = 1e-2

mass = 0.041
inertia_z = 27.8e-6
l_f = 0.04 # wheel base length
l_r = 0.07

# tire model constants
B_f = 3.38
C_f = 1.2
D_f = 0.192

B_r = 4.0
C_r = 1.9
D_r = 0.175

# electric drivetrain and friction constants
C_r0 = 0.0518
C_r2 = 0.00035
C_m1 = 0.287 * 2
C_m2 = 0.054

init_track_constraints = []
track_constraints = []

def transform_track(car_pos):
	phi = car_pos["phi"]
	for i in range(len(init_track_constraints)):
		for element in init_track_constraints[i]:
			diffx = init_track_constraints[i][element]["x"] - car_pos["x"]
			diffy = init_track_constraints[i][element]["y"] - car_pos["y"]
			track_constraints[i][element] = {"x":diffx * cos(phi) + diffy * sin(phi), "y":diffy * cos(phi) - diffx * sin(phi)}


def append_track(left_point, center_point, right_point):
	track_constraints.append({})
	init_track_constraints.append({"left" : left_point, "center" : center_point, "right" : right_point})
def retrieve_track_segment(car_pos, number_of_points, backwards_points, mult):
	result = []
	closest_dist = 99999999
	closest_idx = 0
	for i in range(len(track_constraints)):
		track_element = track_constraints[i]
		dist = sqrt((car_pos["x"] - track_element["center"]["x"])**2 + 1*(car_pos["y"] - track_element["center"]["y"])**2)
		if dist < closest_dist:
			closest_dist = dist
			closest_idx = i
	for i in range(number_of_points):
		idx = (closest_idx + i*mult - backwards_points*mult) % len(track_constraints)
		result.append(track_constraints[idx])

	return result
def list_of_dict_to_list(list_of_dict, key):
	result = []
	for dict_item in list_of_dict:
		result.append(dict_item[key])
	return result

def fit_track_poly(elements, degree):
	x = []
	y = []
	for i in elements:
		x.append(i["x"])
		y.append(i["y"])
	fit = np.polynomial.polynomial.Polynomial.fit(x, y, degree)
	fit = fit.convert(domain=[-1, 1])
	return fit

def generate_track_polys(car_pos, number_of_points, backwards_points, mult):
	track_elements = retrieve_track_segment(car_pos, number_of_points, backwards_points, mult)
	left_poly = fit_track_poly(list_of_dict_to_list(track_elements, 'left'), 4)
	right_poly = fit_track_poly(list_of_dict_to_list(track_elements, 'right'), 4)

	return left_poly, right_poly, track_elements
	

# TEST DATA GENERATE
for track_idx in np.linspace(0, math.pi * 2, 400):
	i = track_idx
	z0 = 1
	if i>math.radians(270):
		z0 = 2
	elif i>math.radians(90):
		z0 = -4

	c0 = 20 + sin(i*4) * z0
	c1 = 18 + sin(i*4) * z0
	c2 = 16 + sin(i*4) * z0
	append_track(
		{"x" : c0 * sin(i), "y" : c0 * cos(i)},
		{"x" : c1 * sin(i), "y" : c1 * cos(i)},
		{"x" : c2 * sin(i), "y" : c2 * cos(i)})

#track_scaling = 5
#for track_idx in range(len(track_mat['trackMobil']['inner'][0][0][0])):
#	outer_x = track_mat['trackMobil']['inner'][0][0][0][track_idx] * track_scaling
#	outer_y = track_mat['trackMobil']['inner'][0][0][1][track_idx] * track_scaling
#	inner_x = track_mat['trackMobil']['outer'][0][0][0][track_idx] * track_scaling
#	inner_y = track_mat['trackMobil']['outer'][0][0][1][track_idx] * track_scaling
#	center_x = track_mat['trackMobil']['center'][0][0][0][track_idx] * track_scaling
#	center_y = track_mat['trackMobil']['center'][0][0][1][track_idx] * track_scaling
#
#	append_track(
#		{"x" : outer_x, "y" : outer_y},
#		{"x" : center_x, "y" : center_y},
#		{"x" : inner_x, "y" : inner_y})

mpc_time = 1.2
mpc_samples = 12

m = ConcreteModel()
m.t = ContinuousSet(bounds=(0, mpc_time))

m.D = Var(m.t) # throttle
m.delta = Var(m.t) # steering angle
m.D_dot = DerivativeVar(m.D)
m.delta_dot = DerivativeVar(m.delta)

m.constr_D_dot_upper = Constraint(m.t, rule=lambda m, t: (m.D_dot[t]) <= 3)
m.constr_D_dot_lower = Constraint(m.t, rule=lambda m, t: (m.D_dot[t]) >= -3)
m.constr_delta_dot_upper = Constraint(m.t, rule=lambda m, t: m.delta_dot[t] <= math.radians(180))
m.constr_delta_dot_lower = Constraint(m.t, rule=lambda m, t: m.delta_dot[t] >= -math.radians(180))

m.obstacle_penalty = Var(m.t)
m.obstacle_caution_penalty = Var(m.t, initialize = 0)

test_Var = 0
world_space_x = 0
world_space_y = 18
world_space_phi = 0

m.x = Var(m.t, initialize=0.0)
m.y = Var(m.t, initialize=0.0)
m.phi = Var(m.t) # heading
m.v_x = Var(m.t)
m.v_y = Var(m.t)
m.omega = DerivativeVar(m.phi) # angular vel

m.x_dot = DerivativeVar(m.x)
m.y_dot = DerivativeVar(m.y)
m.v_x_dot = DerivativeVar(m.v_x)
m.v_y_dot = DerivativeVar(m.v_y)
m.omega_dot = DerivativeVar(m.omega)

m.F_y_f = Var(m.t)
m.F_y_r = Var(m.t)
m.F_x_r = Var(m.t)
m.alpha_f = Var(m.t)
m.alpha_r = Var(m.t)

m.left_track_poly_c0 = Param(mutable=True, default = 0)
m.left_track_poly_c1 = Param(mutable=True, default = 0)
m.left_track_poly_c2 = Param(mutable=True, default = 0)
m.left_track_poly_c3 = Param(mutable=True, default = 0)
m.left_track_poly_c4 = Param(mutable=True, default = 0)
m.right_track_poly_c0 = Param(mutable=True, default = 0)
m.right_track_poly_c1 = Param(mutable=True, default = 0)
m.right_track_poly_c2 = Param(mutable=True, default = 0)
m.right_track_poly_c3 = Param(mutable=True, default = 0)
m.right_track_poly_c4 = Param(mutable=True, default = 0)

m.left_track_constr = Constraint(m.t, rule=lambda m, t: m.y[t] - m.obstacle_penalty[t] <= m.left_track_poly_c0 + m.left_track_poly_c1*m.x[t] + m.left_track_poly_c2*m.x[t]*m.x[t] + m.left_track_poly_c3*m.x[t]*m.x[t]*m.x[t] + m.left_track_poly_c4*m.x[t]*m.x[t]*m.x[t]*m.x[t] - m.obstacle_caution_penalty[t])  
m.right_track_constr = Constraint(m.t, rule=lambda m, t: m.y[t] + m.obstacle_penalty[t] >= m.right_track_poly_c0 + m.right_track_poly_c1*m.x[t] + m.right_track_poly_c2*m.x[t]*m.x[t] + m.right_track_poly_c3*m.x[t]*m.x[t]*m.x[t] + m.right_track_poly_c4*m.x[t]*m.x[t]*m.x[t]*m.x[t] + m.obstacle_caution_penalty[t])

m.x_target = Param(mutable=True, default = 0)
m.y_target = Param(mutable=True, default = 0)

m.D_init = Param(mutable=True, default = 0)
m.delta_init = Param(mutable=True, default = 0)
m.x_init = Param(mutable=True, default = 0)
m.y_init = Param(mutable=True, default = 0)
m.phi_init = Param(mutable=True, default = 0)
m.v_x_init = Param(mutable=True, default = 0)
m.v_y_init = Param(mutable=True, default = 0)
m.omega_init = Param(mutable=True, default = 0)

m.pc = ConstraintList()
m.pc.add(m.D[0]==m.D_init)
m.pc.add(m.delta[0]==m.delta_init)
m.pc.add(m.x[0]==m.x_init)
m.pc.add(m.y[0]==m.y_init)
m.pc.add(m.phi[0]==m.phi_init)
m.pc.add(m.v_x[0]==m.v_x_init)
m.pc.add(m.v_y[0]==m.v_y_init)
m.pc.add(m.omega[0]==m.omega_init)


m.constr_0 = Constraint(m.t, rule=lambda m, t: m.obstacle_caution_penalty[t] >= 0)
m.constr_1 = Constraint(m.t, rule=lambda m, t: m.v_x[t] >= -0.1)
m.constr_2 = Constraint(m.t, rule=lambda m, t: m.D[t] <= 1.0)
m.constr_3 = Constraint(m.t, rule=lambda m, t: m.D[t] >= -2)
m.constr_4 = Constraint(m.t, rule=lambda m, t: m.delta[t] <= math.radians(90))
m.constr_5 = Constraint(m.t, rule=lambda m, t: m.delta[t] >= -math.radians(90))

# kinematic constraints
m.x_dot_ode = Constraint(m.t, rule=lambda m, t: 
	m.x_dot[t] == m.v_x[t] * cos(m.phi[t]) - m.v_y[t] * sin(m.phi[t])
)
m.y_dot_ode = Constraint(m.t, rule=lambda m, t: 
	m.y_dot[t] == m.v_x[t] * sin(m.phi[t]) + m.v_y[t] * cos(m.phi[t])
)
# dynamics constraints
m.v_x_dot_ode = Constraint(m.t, rule=lambda m, t: 
	m.v_x_dot[t] == 1/mass * (m.F_x_r[t] - m.F_y_f[t]*sin(m.delta[t]) + mass*m.v_y[t]*m.omega[t])
)
m.v_y_dot_ode = Constraint(m.t, rule=lambda m, t: 
	m.v_y_dot[t] == 1/mass * (m.F_y_r[t] + m.F_y_f[t]*cos(m.delta[t]) - mass*m.v_x[t]*m.omega[t])
)
m.omega_dot_ode = Constraint(m.t, rule=lambda m, t: 
	m.omega_dot[t] == 1/inertia_z * (m.F_y_f[t]*l_f*cos(m.delta[t]) - m.F_y_r[t]*l_r)
)
# tire / drivetrain dynamics constraints
m.F_y_f_ode = Constraint(m.t, rule=lambda m, t: 
	m.F_y_f[t] == D_f*sin(C_f*atan(B_f*m.alpha_f[t]))
)
m.F_y_r_ode = Constraint(m.t, rule=lambda m, t: 
	m.F_y_r[t] == D_r*sin(C_r*atan(B_r*m.alpha_r[t]))
)
m.F_x_r_ode = Constraint(m.t, rule=lambda m, t: 
	m.F_x_r[t] == (C_m1 - C_m2*m.v_x[t])*m.D[t] - C_r0 - C_r2*pow(m.v_x[t], 2)
)
m.alpha_f_ode = Constraint(m.t, rule=lambda m, t: 
	m.alpha_f[t] == -atan((m.omega[t]*l_f + m.v_y[t]) / (m.v_x[t] + epsilon)) + m.delta[t]
)
m.alpha_r_ode = Constraint(m.t, rule=lambda m, t: 
	m.alpha_r[t] == atan((m.omega[t]*l_r - m.v_y[t]) / (m.v_x[t] + epsilon))
)

#TransformationFactory('dae.finite_difference').apply_to(m, wrt=m.t, nfe=mpc_samples)
TransformationFactory('dae.collocation').apply_to(m, wrt=m.t, nfe=mpc_samples, ncp=1, scheme='LAGRANGE-RADAU')

def cost_function_integral(m, t):
	return 0.25*(m.x[t] - m.x_target)**2 + 0.25*(m.y[t] - m.y_target)**2 + 900*(m.obstacle_penalty[t])**2 + 200*(m.obstacle_caution_penalty[t]-1)**2 + -10*(m.v_x[t])**2 #+ 1*(m.delta[t])**2

m.integral = Integral(m.t, wrt=m.t, rule=cost_function_integral)
m.obj = Objective(expr=m.integral)

fig, (ax1, ax2, ax3) = plt.subplots(3,1)
#fig, ax1 = plt.subplots(1,1)

# xy, left_poly, right_poly, left_track, right_track, center, actual_left_Track, actual_right_track
lines = [(ax1.plot([], [], linestyle='--', marker='x', color='g')[0]), (ax1.plot([], [], lw=2)[0]), (ax1.plot([], [], lw=2)[0]), (ax1.plot([], [], '--bo')[0]), (ax1.plot([], [], '--bo')[0]), (ax1.plot([], [], linestyle='None', marker='x', color='pink')[0]), (ax1.plot([], [], linestyle='--', color = 'silver')[0]), (ax1.plot([], [], linestyle='--', color = 'silver')[0]), 
(ax2.plot([], [], linestyle='-', color = 'black')[0]), (ax2.plot([], [], linestyle=':', color = 'red')[0]), (ax2.plot([], [], linestyle=':', color = 'blue')[0]), (ax2.plot([], [], linestyle='-', color = 'limegreen')[0]), (ax1.plot([], [], linestyle='--', color = 'red')[0])]
#ax1.plot([0,10], [0,0], linestyle='--', color = 'black')

lines_len_init =len(lines)-1

plot_list = {}
plot_params = ['delta', 'D', 'v_x', 'v_y']

axis_lims_inner = 8
axis_lims_outer = 22
ax1.set_xlim(-axis_lims_inner, axis_lims_inner)
ax1.set_ylim(-axis_lims_inner, axis_lims_inner)
ax1.set_aspect('equal', 'box')
ax2.set_xlim(-axis_lims_outer, axis_lims_outer)
ax2.set_ylim(-axis_lims_outer, axis_lims_outer)
ax2.set_aspect('equal', 'box')
ax3.set_xlim(0, mpc_time)
ax3.set_ylim(-2, 2)
ax3.set_aspect(0.1)

world_space_x_list = []
world_space_y_list = []

def solve(solve_idx):
	global world_space_x, world_space_y, world_space_phi, car_arrow
	transform_track({"x": world_space_x, "y": world_space_y, "phi": world_space_phi})

	left_track_poly, right_track_poly, track_elements = generate_track_polys({"x":m.x[0](), "y":m.y[0]()}, 6, 1, 8)
	target_track_element = track_elements[-1]["center"]
	m.x_target = target_track_element["x"]
	m.y_target = target_track_element["y"]

	m.left_track_poly_c0 = left_track_poly.coef[0]
	m.left_track_poly_c1 = left_track_poly.coef[1]
	m.left_track_poly_c2 = left_track_poly.coef[2]
	m.left_track_poly_c3 = left_track_poly.coef[3]
	m.left_track_poly_c4 = left_track_poly.coef[4]

	m.right_track_poly_c0 = right_track_poly.coef[0]
	m.right_track_poly_c1 = right_track_poly.coef[1]
	m.right_track_poly_c2 = right_track_poly.coef[2]
	m.right_track_poly_c3 = right_track_poly.coef[3]
	m.right_track_poly_c4 = right_track_poly.coef[4]

	track_x_space = np.linspace(m.x[0]()-15,m.x[0]()+15,300)
	lines[1].set_data(track_x_space, left_track_poly(track_x_space))
	lines[2].set_data(track_x_space, right_track_poly(track_x_space))
	lines[3].set_data(list_of_dict_to_list(list_of_dict_to_list(track_elements, 'left'), 'x'), list_of_dict_to_list(list_of_dict_to_list(track_elements, 'left'), 'y'))
	lines[4].set_data(list_of_dict_to_list(list_of_dict_to_list(track_elements, 'right'), 'x'), list_of_dict_to_list(list_of_dict_to_list(track_elements, 'right'), 'y'))
	lines[5].set_data(list_of_dict_to_list(list_of_dict_to_list(track_elements, 'center'), 'x'), list_of_dict_to_list(list_of_dict_to_list(track_elements, 'center'), 'y'))
	lines[6].set_data(list_of_dict_to_list(list_of_dict_to_list(track_constraints, 'left'), 'x'), list_of_dict_to_list(list_of_dict_to_list(track_constraints, 'left'), 'y'))
	lines[7].set_data(list_of_dict_to_list(list_of_dict_to_list(track_constraints, 'right'), 'x'), list_of_dict_to_list(list_of_dict_to_list(track_constraints, 'right'), 'y'))
	lines[8].set_data(world_space_x_list, world_space_y_list)
	lines[9].set_data(list_of_dict_to_list(list_of_dict_to_list(init_track_constraints, 'left'), 'x'), list_of_dict_to_list(list_of_dict_to_list(init_track_constraints, 'left'), 'y'))
	lines[10].set_data(list_of_dict_to_list(list_of_dict_to_list(init_track_constraints, 'right'), 'x'), list_of_dict_to_list(list_of_dict_to_list(init_track_constraints, 'right'), 'y'))
	lines[11].set_data([[world_space_x + -1*cos(world_space_phi), world_space_x + cos(world_space_phi)], [world_space_y + -1*sin(world_space_phi), world_space_y + sin(world_space_phi)]])
	lines[12].set_data([m.v_x_init(),0], [m.v_y_init(),0])

	if len(world_space_x_list) > 100:
		for i in range(10):
			world_space_x_list.pop(0)
			world_space_y_list.pop(0)

	#print(m.x_target(), m.y_target())
	#caution_sum = 0
	#for t in m.t:
	#	caution_sum = caution_sum + m.obstacle_caution_penalty[t]()

	last_solve_time = time.time()
	solver.solve(m, tee=False)
	print(1/(time.time() - last_solve_time))

	init_next_idx = m.t[2]

	world_space_phi = m.phi[init_next_idx]() + world_space_phi
	world_space_x = world_space_x + (m.x[init_next_idx]() * cos(world_space_phi) + m.y[init_next_idx]() * sin(world_space_phi))
	world_space_y = world_space_y + (m.y[init_next_idx]() * cos(world_space_phi) + m.x[init_next_idx]() * sin(world_space_phi))
	world_space_x_list.append(world_space_x)
	world_space_y_list.append(world_space_y)

	m.D_init = m.D[init_next_idx]()
	m.delta_init = m.delta[init_next_idx]()
	#m.x_init = m.x[init_next_idx]()
	#m.y_init = m.y[init_next_idx]()
	#m.phi_init = m.phi[init_next_idx]()
	m.v_x_init = m.v_x[init_next_idx]()
	m.v_y_init = m.v_y[init_next_idx]()
	m.omega_init = m.omega[init_next_idx]()

	t_array = np.array([t for t in m.t]).tolist()
	x_array = np.array([m.x[t]() for t in m.t]).tolist()
	y_array = np.array([m.y[t]() for t in m.t]).tolist()

	lines[0].set_data(x_array, y_array)

	for element in plot_params:
		uniq_mul = 1
		if element == 'v_x':
			uniq_mul = 0.2
		elif element == 'delta':
			uniq_mul = math.degrees(0.5)
		plot_list[element] = np.array([getattr(m, element)[t]()*uniq_mul for t in m.t]).tolist()
	line_idx = lines_len_init
	for key in plot_list:
		if len(lines)-1 <= line_idx:
			new_line, = ax3.plot([], [], label=key, linestyle='--')
			lines.append(new_line)
		line_idx = line_idx + 1
		lines[line_idx].set_data(t_array,plot_list[key])

	plt.legend(loc='right')

	return lines

anim = animation.FuncAnimation(fig, solve, 100000, interval=1, blit=True, repeat=False)
plt.show()