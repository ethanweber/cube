import numpy as np
from pydrake.math import sin, cos
from pydrake.all import *

from constraints import *
from dynamics_nd import *

# dimension 3 works well
# initial states in the code
def swing_up(dim=2):

    x = 0.; y = 0.0; theta = 0.;
    initial_state = get_nd_state((x,y,theta,0,0,0,0,0), dim)
    x = -0.5; y = .5*(2**.5)-.5; theta = np.pi/4.0;
    final_state = get_nd_state((x,y,theta,0,0,0,0,0), dim)
    min_time = 0.5; max_time = 15.; max_torque = 1000.0

    print("Initial State: {}".format(initial_state))
    print("Final State: {}".format(final_state))

    # a few checks
    assert(len(initial_state) == len(final_state))
    assert(min_time <= max_time)

    # some values that can be changed if desired
    N = 50 # number knot points
    dynamics_error_thresh = 0.01 # error thresh on xdot = f(x,u,f)
    floor_offset = -.01 # used to allow a little penitration
    final_state_error_thresh = 0.01 # final state error thresh
    max_ground_force = 100
    # impose contraint to stay on the ground
    stay_on_ground = True
    stay_on_ground_tolerance = 0.1 # tolerance for leaving the ground if contstraint is used
    complimentarity_constraint_thresh = 0.0001

    fix_corner_on_ground = True
    corner_fix = [0, -0.5] # corner index, location

    mu = 0.1 # friction force

    # use SNOPT in Drake
    mp = MathematicalProgram()

    # state length
    state_len = len(initial_state)

    # total time used (assuming equal time steps)
    time_used = mp.NewContinuousVariables(1, "time_used")
    dt = time_used/(N+1)

    # input torque decision variables
    u = mp.NewContinuousVariables(1, "u_%d" % 0) # only one input for the cube
    u_over_time = u
    for k in range(1,N):
        u = mp.NewContinuousVariables(1, "u_%d" % k)
        u_over_time = np.vstack((u_over_time, u))
    total_u = u_over_time

    # contact force decision variables
    f = mp.NewContinuousVariables(8, "f_%d" % 0) # only one input for the cube
    f_over_time = f
    for k in range(1,N):
        f = mp.NewContinuousVariables(8, "f_%d" % k)
        f_over_time = np.vstack((f_over_time, f))
    total_f = f_over_time

    # state decision variables
    x = mp.NewContinuousVariables(state_len, "x_%d" % 0) # for both input thrusters
    x_over_time = x
    for k in range(1,N+1):
        x = mp.NewContinuousVariables(state_len, "x_%d" % k)
        x_over_time = np.vstack((x_over_time, x))
    total_x = x_over_time

    # impose dynamic constraints
    for n in range(N):
        state_next = total_x[n+1]
        dynamic_state_next = total_x[n,:] + get_nd_dynamics(total_x[n,:], total_u[n,:], total_f[n,:], dim)*dt
        # make sure the actual and predicted align to follow dynamics
        for j in range(state_len):
            state_error = state_next[j] - dynamic_state_next[j]
            mp.AddConstraint(state_error <= dynamics_error_thresh)
            mp.AddConstraint(state_error >= -dynamics_error_thresh)

    # can't penitrate the floor and can't leave the floor
    for n in range(N):
        distances = get_corner_distances(total_x[n,:], dim)

        mp.AddConstraint(distances[0] >= floor_offset)
        mp.AddConstraint(distances[1] >= floor_offset)
        mp.AddConstraint(distances[2] >= floor_offset)
        mp.AddConstraint(distances[3] >= floor_offset)

        # don't leave the ground if specified
        if stay_on_ground == True:
            mp.AddConstraint(distances[0] <= np.sqrt(2)+stay_on_ground_tolerance)
            mp.AddConstraint(distances[1] <= np.sqrt(2)+stay_on_ground_tolerance)
            mp.AddConstraint(distances[2] <= np.sqrt(2)+stay_on_ground_tolerance)
            mp.AddConstraint(distances[3] <= np.sqrt(2)+stay_on_ground_tolerance)

        if fix_corner_on_ground == True:
            x_pos = get_corner_x_positions(total_x[n,:], dim)
            # make left corner on the ground
            mp.AddConstraint(distances[0] == 0.0)
            num, loc = corner_fix
            mp.AddConstraint(x_pos[num] == loc)

    # ground forces can't pull on the ground
    for n in range(N):
        force = total_f[n]
        for j in range(8):
            mp.AddConstraint(force[j] <= max_ground_force)
            mp.AddConstraint(force[j] >= 0)

    # add complimentary constraint
    for n in range(N):
        force = total_f[n]
        state = total_x[n]
        theta = state[dim]

        distances = get_corner_distances(total_x[n+1,:], dim)

        s = sin(theta)
        c = cos(theta)

        z_0 = force[0]*c + force[1]*s
        z_1 = - force[2]*s + force[3]*c
        z_2 = - force[4]*c - force[5]*s
        z_3 = force[6]*s - force[7]*c

        xy_0 = - force[0]*s + force[1]*c
        xy_1 = - force[2]*c - force[3]*s
        xy_2 = force[4]*s - force[5]*c
        xy_3 = force[6]*c + force[7]*s

        mp.AddConstraint(xy_0 <= z_0*mu)
        mp.AddConstraint(xy_0 >= -z_0*mu)
        mp.AddConstraint(xy_1 <= z_1*mu)
        mp.AddConstraint(xy_1 >= -z_1*mu)
        mp.AddConstraint(xy_2 <= z_2*mu)
        mp.AddConstraint(xy_2 >= -z_2*mu)
        mp.AddConstraint(xy_3 <= z_3*mu)
        mp.AddConstraint(xy_3 >= -z_3*mu)

        # vector_0 = force[0] * force[1]
        # vector_1 = force[2] * force[3]
        # vector_2 = force[4] * force[5]
        # vector_3 = force[6] * force[7]
        # val = np.asarray([vector_0, vector_1, vector_2, vector_3])
        # mp.AddConstraint(val.dot(distances) <= complimentarity_constraint_thresh)
        # mp.AddConstraint(val.dot(distances) >= -complimentarity_constraint_thresh)

        val_0 = np.asarray([force[0], force[2], force[4], force[6]])
        val_1 = np.asarray([force[1], force[3], force[5], force[7]])
        mp.AddConstraint(val_0.dot(distances) <= complimentarity_constraint_thresh)
        mp.AddConstraint(val_0.dot(distances) >= -complimentarity_constraint_thresh)
        mp.AddConstraint(val_1.dot(distances) <= complimentarity_constraint_thresh)
        mp.AddConstraint(val_1.dot(distances) >= -complimentarity_constraint_thresh)

    # initial state, no state error allowed
    for i in range(state_len):
        initial_state_error = x_over_time[0,i] - initial_state[i]
        mp.AddConstraint(initial_state_error == 0.0)

    # don't care about final wheel angle, so skip restriction in it
    state_indices = [i for i in range(0, state_len)]
    a = state_indices[0:state_len/2-1] + state_indices[state_len/2:]
    for i in a:
        # final
        final_state_error = x_over_time[-1,i] - final_state[i]
        mp.AddConstraint(final_state_error <= final_state_error_thresh)
        mp.AddConstraint(final_state_error >= -final_state_error_thresh)

    # add time constraint
    mp.AddConstraint(time_used[0] >= min_time)
    mp.AddConstraint(time_used[0] <= max_time)

    # add torque constraints
    for n in range(N):
        mp.AddConstraint(u_over_time[n,0] <= max_torque)
        mp.AddConstraint(u_over_time[n,0] >= -max_torque)

    # try to keep the velocity of the wheel in the correct direction
    # mp.AddLinearCost(x_over_time[:,-1].sum())

    # mp.AddLinearCost(-x_over_time[:,1].sum())
    # mp.AddLinearCost(-x_over_time[N//2,1])

    # minimize input
    # mp.AddQuadraticCost(u_over_time[:,0].dot(u_over_time[:,0]))

    # minimize the time used
    # mp.AddLinearCost(time_used[0])

    print "Number of decision vars", mp.num_vars()
    print(mp.Solve())

    trajectory = mp.GetSolution(x_over_time)
    input_trajectory = mp.GetSolution(u_over_time)
    force_trajectory = mp.GetSolution(f_over_time)
    t = mp.GetSolution(time_used)
    time_array = np.arange(0.0, t, t/(N+1))

    return trajectory, input_trajectory, force_trajectory, time_array

# periodic motion
def periodic_motion(dim=2):

    # print("Initial State: {}".format(initial_state))
    # print("Final State: {}".format(final_state))

    # a few checks
    min_time = 0.5
    max_time = 15.0
    max_torque = 1000.0

    # some values that can be changed if desired
    N = 50 # number knot points
    dynamics_error_thresh = 0.01 # error thresh on xdot = f(x,u,f)
    floor_offset = -.01 # used to allow a little penitration
    final_state_error_thresh = 0.01 # final state error thresh
    max_ground_force = 100
    # impose contraint to stay on the ground
    stay_on_ground = True
    stay_on_ground_tolerance = 0.1 # tolerance for leaving the ground if contstraint is used
    complimentarity_constraint_thresh = 0.0001

    fix_corner_on_ground = True
    corner_fix = [0, -0.5] # corner index, location

    mu = 0.1 # friction force

    # use SNOPT in Drake
    mp = MathematicalProgram()

    # state length
    state_len = 8

    # total time used (assuming equal time steps)
    time_used = mp.NewContinuousVariables(1, "time_used")
    dt = time_used/(N+1)

    # input torque decision variables
    u = mp.NewContinuousVariables(1, "u_%d" % 0) # only one input for the cube
    u_over_time = u
    for k in range(1,N):
        u = mp.NewContinuousVariables(1, "u_%d" % k)
        u_over_time = np.vstack((u_over_time, u))
    total_u = u_over_time

    # contact force decision variables
    f = mp.NewContinuousVariables(8, "f_%d" % 0) # only one input for the cube
    f_over_time = f
    for k in range(1,N):
        f = mp.NewContinuousVariables(8, "f_%d" % k)
        f_over_time = np.vstack((f_over_time, f))
    total_f = f_over_time

    # state decision variables
    x = mp.NewContinuousVariables(state_len, "x_%d" % 0) # for both input thrusters
    x_over_time = x
    for k in range(1,N+1):
        x = mp.NewContinuousVariables(state_len, "x_%d" % k)
        x_over_time = np.vstack((x_over_time, x))
    total_x = x_over_time

    # impose dynamic constraints
    for n in range(N):
        state_next = total_x[n+1]
        dynamic_state_next = total_x[n,:] + get_nd_dynamics(total_x[n,:], total_u[n,:], total_f[n,:], dim)*dt
        # make sure the actual and predicted align to follow dynamics
        for j in range(state_len):
            state_error = state_next[j] - dynamic_state_next[j]
            mp.AddConstraint(state_error <= dynamics_error_thresh)
            mp.AddConstraint(state_error >= -dynamics_error_thresh)

    # can't penitrate the floor and can't leave the floor
    for n in range(N):
        distances = get_corner_distances(total_x[n,:], dim)

        mp.AddConstraint(distances[0] >= floor_offset)
        mp.AddConstraint(distances[1] >= floor_offset)
        mp.AddConstraint(distances[2] >= floor_offset)
        mp.AddConstraint(distances[3] >= floor_offset)

        # don't leave the ground if specified
        if stay_on_ground == True:
            mp.AddConstraint(distances[0] <= np.sqrt(2)+stay_on_ground_tolerance)
            mp.AddConstraint(distances[1] <= np.sqrt(2)+stay_on_ground_tolerance)
            mp.AddConstraint(distances[2] <= np.sqrt(2)+stay_on_ground_tolerance)
            mp.AddConstraint(distances[3] <= np.sqrt(2)+stay_on_ground_tolerance)

        if fix_corner_on_ground == True:
            x_pos = get_corner_x_positions(total_x[n,:], dim)
            # make left corner on the ground
            mp.AddConstraint(distances[0] == 0.0)
            num, loc = corner_fix
            mp.AddConstraint(x_pos[num] == loc)

    # ground forces can't pull on the ground
    for n in range(N):
        force = total_f[n]
        for j in range(8):
            mp.AddConstraint(force[j] <= max_ground_force)
            mp.AddConstraint(force[j] >= 0)

    # add complimentary constraint
    for n in range(N):
        force = total_f[n]
        state = total_x[n]
        theta = state[dim]

        distances = get_corner_distances(total_x[n+1,:], dim)

        s = sin(theta)
        c = cos(theta)

        z_0 = force[0]*c + force[1]*s
        z_1 = - force[2]*s + force[3]*c
        z_2 = - force[4]*c - force[5]*s
        z_3 = force[6]*s - force[7]*c

        xy_0 = - force[0]*s + force[1]*c
        xy_1 = - force[2]*c - force[3]*s
        xy_2 = force[4]*s - force[5]*c
        xy_3 = force[6]*c + force[7]*s

        mp.AddConstraint(xy_0 <= z_0*mu)
        mp.AddConstraint(xy_0 >= -z_0*mu)
        mp.AddConstraint(xy_1 <= z_1*mu)
        mp.AddConstraint(xy_1 >= -z_1*mu)
        mp.AddConstraint(xy_2 <= z_2*mu)
        mp.AddConstraint(xy_2 >= -z_2*mu)
        mp.AddConstraint(xy_3 <= z_3*mu)
        mp.AddConstraint(xy_3 >= -z_3*mu)

        val_0 = np.asarray([force[0], force[2], force[4], force[6]])
        val_1 = np.asarray([force[1], force[3], force[5], force[7]])
        mp.AddConstraint(val_0.dot(distances) <= complimentarity_constraint_thresh)
        mp.AddConstraint(val_0.dot(distances) >= -complimentarity_constraint_thresh)
        mp.AddConstraint(val_1.dot(distances) <= complimentarity_constraint_thresh)
        mp.AddConstraint(val_1.dot(distances) >= -complimentarity_constraint_thresh)

    # initial
    mp.AddConstraint(x_over_time[0,0] == 0.0) # x
    mp.AddConstraint(x_over_time[0,1] == 0.0) # y
    mp.AddConstraint(x_over_time[0,2] == 0.0) # 0 angle
    mp.AddConstraint(x_over_time[0,3] == 0.0) # alpha angle
    # final
    mp.AddConstraint(x_over_time[-1,0] == -1.0) # x
    mp.AddConstraint(x_over_time[-1,1] == 0.0) # y
    mp.AddConstraint(x_over_time[-1,2] == np.pi/2.0) # y

    # connection constraint, don't care about inner wheel angle
    for i in [4,5,6,7]:
        mp.AddConstraint(x_over_time[0,i] == x_over_time[-1,i])

    # add time constraint
    mp.AddConstraint(time_used[0] >= min_time)
    mp.AddConstraint(time_used[0] <= max_time)

    # add torque constraints
    for n in range(N):
        mp.AddConstraint(u_over_time[n,0] <= max_torque)
        mp.AddConstraint(u_over_time[n,0] >= -max_torque)

    # minimize input
    # mp.AddQuadraticCost(u_over_time[:,0].dot(u_over_time[:,0]))

    # maximize velocity in correct direction (left)
    mp.AddLinearCost(-x_over_time[:,4].sum())

    # minimize the time
    # mp.AddLinearCost(time_used[0])

    print "Number of decision vars", mp.num_vars()
    print(mp.Solve())

    trajectory = mp.GetSolution(x_over_time)
    input_trajectory = mp.GetSolution(u_over_time)
    force_trajectory = mp.GetSolution(f_over_time)
    t = mp.GetSolution(time_used)
    time_array = np.arange(0.0, t, t/(N+1))

    return trajectory, input_trajectory, force_trajectory, time_array

# upright controller
def qp_controller(current_state, desired_state, dt, dim=2):
    """This is the controller that returns the torque for the desired state.

    Keyword arguments:
    current_state -- current state
    desired_state -- desired state
    dt -- timestep
    dim -- state dimension (default 2)
    """

    # torque PD controller values
    wheel_kp = 50.0
    wheel_kd = 10.0
    max_torque = 20.0

    # cost on obtaining next state and velocity
    kp = 0.0
    kd = 1.0

    # half state length
    hl = len(current_state) / 2

    mp = MathematicalProgram()

    x = mp.NewContinuousVariables(len(current_state), "x")
    u = mp.NewContinuousVariables(1, "u")
    force = mp.NewContinuousVariables(8, "force")

    # set the initial state
    set_initial_state(mp, x, current_state, dim)
    # enforce the dynamics with linearized theta
    state = x + get_nd_dynamics(x, u, force, dim, current_state[dim])*dt

    # stay on floor
    # add_floor_constraint(mp, state, dim)
    # for corner to ground
    # fix_corner_to_ground(mp, state, 0, -0.5, dim)
    # don't pull on ground
    dont_pull_on_ground(mp, force, dim)
    # bounded to not leave the ground
    # stay_on_ground(mp, state, dim)
    # only force when on ground
    complimentarity_constraint(mp, state, force, dim)

    # linearize theta to set this cost
    add_corner_cost(mp, state, 0, -0.5, dim, current_state[dim])

    # unpack the states
    x_s = state[0]
    y = state[1]
    theta = state[dim]
    alpha = state[hl-1]
    xdot = state[0+hl]
    ydot = state[1+hl]
    theta_dot = state[dim+hl]
    alpha_dot = state[-1]

    # unpack the desired states
    x_des = desired_state[0]
    y_des = desired_state[1]
    theta_des = desired_state[dim]
    alpha_des = desired_state[hl-1]
    xdot_des = desired_state[0+hl]
    ydot_des = desired_state[1+hl]
    theta_dot_des = desired_state[dim+hl]
    alpha_dot_des = desired_state[-1]

    # current_pos = np.asarray([x_s,y,theta,alpha])
    # des_pos = np.asarray([x_des,y_des,theta_des,alpha_des])
    # pos_diff = current_pos - des_pos
    current_pos = np.asarray([x_s,y,theta,0])
    des_pos = np.asarray([x_des,y_des,theta_des,0])
    pos_diff = current_pos - des_pos

    # current_vel = np.asarray([xdot,ydot,theta_dot,alpha_dot])
    # des_vel = np.asarray([xdot_des,ydot_des,theta_dot_des,alpha_dot_des])
    # vel_diff = current_vel - des_vel
    current_vel = np.asarray([xdot,ydot,theta_dot,0])
    des_vel = np.asarray([xdot_des,ydot_des,theta_dot_des,0])
    vel_diff = current_vel - des_vel

    pos = pos_diff.dot(pos_diff)
    vel = vel_diff.dot(vel_diff)

    mp.AddQuadraticCost(kp*pos)
    mp.AddQuadraticCost(kd*vel)

    # torque PD controller
    input_torque = wheel_kp*(current_state[dim] - np.pi/4.0) + wheel_kd*current_state[dim+hl]
    input_torque = np.clip(input_torque, -max_torque, max_torque)
    mp.AddConstraint(u[0] == input_torque)

    sol = mp.Solve()
    # print(sol)

    my_torque = mp.GetSolution(u)
    my_force = mp.GetSolution(force)
    my_start = mp.GetSolution(x)

    return my_start, my_torque, my_force

# general controller
def compute_optimal_control(initial_state, final_state, min_time, max_time, max_torque, dim=2):

    print("Initial State: {}".format(initial_state))
    print("Final State: {}".format(final_state))

    # a few checks
    assert(len(initial_state) == len(final_state))
    assert(min_time <= max_time)

    # some values that can be changed if desired
    N = 50 # number knot points
    dynamics_error_thresh = 0.01 # error thresh on xdot = f(x,u,f)
    floor_offset = -.01 # used to allow a little penitration
    final_state_error_thresh = 0.01 # final state error thresh
    max_ground_force = 100
    # impose contraint to stay on the ground
    stay_on_ground = True
    stay_on_ground_tolerance = 0.1 # tolerance for leaving the ground if contstraint is used
    complimentarity_constraint_thresh = 0.0001

    fix_corner_on_ground = True
    corner_fix = [0, -0.5] # corner index, location

    mu = 0.1 # friction force

    # use SNOPT in Drake
    mp = MathematicalProgram()

    # state length
    state_len = len(initial_state)

    # total time used (assuming equal time steps)
    time_used = mp.NewContinuousVariables(1, "time_used")
    dt = time_used/(N+1)

    # input torque decision variables
    u = mp.NewContinuousVariables(1, "u_%d" % 0) # only one input for the cube
    u_over_time = u
    for k in range(1,N):
        u = mp.NewContinuousVariables(1, "u_%d" % k)
        u_over_time = np.vstack((u_over_time, u))
    total_u = u_over_time

    # contact force decision variables
    f = mp.NewContinuousVariables(8, "f_%d" % 0) # only one input for the cube
    f_over_time = f
    for k in range(1,N):
        f = mp.NewContinuousVariables(8, "f_%d" % k)
        f_over_time = np.vstack((f_over_time, f))
    total_f = f_over_time

    # state decision variables
    x = mp.NewContinuousVariables(state_len, "x_%d" % 0) # for both input thrusters
    x_over_time = x
    for k in range(1,N+1):
        x = mp.NewContinuousVariables(state_len, "x_%d" % k)
        x_over_time = np.vstack((x_over_time, x))
    total_x = x_over_time

    # impose dynamic constraints
    for n in range(N):
        state_next = total_x[n+1]
        dynamic_state_next = total_x[n,:] + get_nd_dynamics(total_x[n,:], total_u[n,:], total_f[n,:], dim)*dt
        # make sure the actual and predicted align to follow dynamics
        for j in range(state_len):
            state_error = state_next[j] - dynamic_state_next[j]
            mp.AddConstraint(state_error <= dynamics_error_thresh)
            mp.AddConstraint(state_error >= -dynamics_error_thresh)

    # can't penitrate the floor and can't leave the floor
    for n in range(N):
        distances = get_corner_distances(total_x[n,:], dim)

        mp.AddConstraint(distances[0] >= floor_offset)
        mp.AddConstraint(distances[1] >= floor_offset)
        mp.AddConstraint(distances[2] >= floor_offset)
        mp.AddConstraint(distances[3] >= floor_offset)

        # don't leave the ground if specified
        if stay_on_ground == True:
            mp.AddConstraint(distances[0] <= np.sqrt(2)+stay_on_ground_tolerance)
            mp.AddConstraint(distances[1] <= np.sqrt(2)+stay_on_ground_tolerance)
            mp.AddConstraint(distances[2] <= np.sqrt(2)+stay_on_ground_tolerance)
            mp.AddConstraint(distances[3] <= np.sqrt(2)+stay_on_ground_tolerance)

        if fix_corner_on_ground == True:
            x_pos = get_corner_x_positions(total_x[n,:], dim)
            # make left corner on the ground
            mp.AddConstraint(distances[0] == 0.0)
            num, loc = corner_fix
            mp.AddConstraint(x_pos[num] == loc)

    # ground forces can't pull on the ground
    for n in range(N):
        force = total_f[n]
        for j in range(8):
            mp.AddConstraint(force[j] <= max_ground_force)
            mp.AddConstraint(force[j] >= 0)

    # add complimentary constraint
    for n in range(N):
        force = total_f[n]
        state = total_x[n]
        theta = state[dim]

        distances = get_corner_distances(total_x[n+1,:], dim)

        s = sin(theta)
        c = cos(theta)

        z_0 = force[0]*c + force[1]*s
        z_1 = - force[2]*s + force[3]*c
        z_2 = - force[4]*c - force[5]*s
        z_3 = force[6]*s - force[7]*c

        xy_0 = - force[0]*s + force[1]*c
        xy_1 = - force[2]*c - force[3]*s
        xy_2 = force[4]*s - force[5]*c
        xy_3 = force[6]*c + force[7]*s

        mp.AddConstraint(xy_0 <= z_0*mu)
        mp.AddConstraint(xy_0 >= -z_0*mu)
        mp.AddConstraint(xy_1 <= z_1*mu)
        mp.AddConstraint(xy_1 >= -z_1*mu)
        mp.AddConstraint(xy_2 <= z_2*mu)
        mp.AddConstraint(xy_2 >= -z_2*mu)
        mp.AddConstraint(xy_3 <= z_3*mu)
        mp.AddConstraint(xy_3 >= -z_3*mu)

        # vector_0 = force[0] * force[1]
        # vector_1 = force[2] * force[3]
        # vector_2 = force[4] * force[5]
        # vector_3 = force[6] * force[7]
        # val = np.asarray([vector_0, vector_1, vector_2, vector_3])
        # mp.AddConstraint(val.dot(distances) <= complimentarity_constraint_thresh)
        # mp.AddConstraint(val.dot(distances) >= -complimentarity_constraint_thresh)

        val_0 = np.asarray([force[0], force[2], force[4], force[6]])
        val_1 = np.asarray([force[1], force[3], force[5], force[7]])
        mp.AddConstraint(val_0.dot(distances) <= complimentarity_constraint_thresh)
        mp.AddConstraint(val_0.dot(distances) >= -complimentarity_constraint_thresh)
        mp.AddConstraint(val_1.dot(distances) <= complimentarity_constraint_thresh)
        mp.AddConstraint(val_1.dot(distances) >= -complimentarity_constraint_thresh)

    # initial state, no state error allowed
    for i in range(state_len):
        initial_state_error = x_over_time[0,i] - initial_state[i]
        mp.AddConstraint(initial_state_error == 0.0)

    # don't care about final wheel angle, so skip restriction in it
    state_indices = [i for i in range(0, state_len)]
    a = state_indices[0:state_len/2-1] + state_indices[state_len/2:]
    for i in a:
        # final
        final_state_error = x_over_time[-1,i] - final_state[i]
        mp.AddConstraint(final_state_error <= final_state_error_thresh)
        mp.AddConstraint(final_state_error >= -final_state_error_thresh)

    # add time constraint
    mp.AddConstraint(time_used[0] >= min_time)
    mp.AddConstraint(time_used[0] <= max_time)

    # add torque constraints
    for n in range(N):
        mp.AddConstraint(u_over_time[n,0] <= max_torque)
        mp.AddConstraint(u_over_time[n,0] >= -max_torque)

    # try to keep the velocity of the wheel in the correct direction
    # mp.AddLinearCost(x_over_time[:,-1].sum())

    # mp.AddLinearCost(-x_over_time[:,1].sum())
    # mp.AddLinearCost(-x_over_time[N//2,1])

    print "Number of decision vars", mp.num_vars()
    print(mp.Solve())

    trajectory = mp.GetSolution(x_over_time)
    input_trajectory = mp.GetSolution(u_over_time)
    force_trajectory = mp.GetSolution(f_over_time)
    t = mp.GetSolution(time_used)
    time_array = np.arange(0.0, t, t/(N+1))

    return trajectory, input_trajectory, force_trajectory, time_array

# TODO make this control work
def mblock_example(initial_state, final_state, min_time, max_time, max_torque, dim=2):

    print("Initial State: {}".format(initial_state))
    print("Final State: {}".format(final_state))

    # a few checks
    assert(len(initial_state) == len(final_state))
    assert(min_time <= max_time)

    # some values that can be changed if desired
    N = 50 # number knot points
    dynamics_error_thresh = 0.01 # error thresh on xdot = f(x,u,f)
    floor_offset = -.01 # used to allow a little penitration
    final_state_error_thresh = 0.01 # final state error thresh
    max_ground_force = 100
    # impose contraint to stay on the ground
    stay_on_ground = False
    stay_on_ground_tolerance = 0.1 # tolerance for leaving the ground if contstraint is used
    complimentarity_constraint_thresh = 0.0001

    fix_corner_on_ground = True
    corner_fix = [2, 0.5] # corner index, location

    mu = 0.1 # friction force

    # use SNOPT in Drake
    mp = MathematicalProgram()

    # state length
    state_len = len(initial_state)

    # total time used (assuming equal time steps)
    time_used = mp.NewContinuousVariables(1, "time_used")
    dt = time_used/(N+1)

    # input torque decision variables
    u = mp.NewContinuousVariables(1, "u_%d" % 0) # only one input for the cube
    u_over_time = u
    for k in range(1,N):
        u = mp.NewContinuousVariables(1, "u_%d" % k)
        u_over_time = np.vstack((u_over_time, u))
    total_u = u_over_time

    # contact force decision variables
    f = mp.NewContinuousVariables(8, "f_%d" % 0) # only one input for the cube
    f_over_time = f
    for k in range(1,N):
        f = mp.NewContinuousVariables(8, "f_%d" % k)
        f_over_time = np.vstack((f_over_time, f))
    total_f = f_over_time

    # state decision variables
    x = mp.NewContinuousVariables(state_len, "x_%d" % 0) # for both input thrusters
    x_over_time = x
    for k in range(1,N+1):
        x = mp.NewContinuousVariables(state_len, "x_%d" % k)
        x_over_time = np.vstack((x_over_time, x))
    total_x = x_over_time

    # impose dynamic constraints
    for n in range(N):
        state_next = total_x[n+1]
        dynamic_state_next = total_x[n,:] + get_nd_dynamics(total_x[n,:], total_u[n,:], total_f[n,:], dim)*dt
        # make sure the actual and predicted align to follow dynamics
        for j in range(state_len):
            state_error = state_next[j] - dynamic_state_next[j]
            mp.AddConstraint(state_error <= dynamics_error_thresh)
            mp.AddConstraint(state_error >= -dynamics_error_thresh)

    # can't penitrate the floor and can't leave the floor
    for n in range(N):
        distances = get_corner_distances(total_x[n,:], dim)

        mp.AddConstraint(distances[0] >= floor_offset)
        mp.AddConstraint(distances[1] >= floor_offset)
        mp.AddConstraint(distances[2] >= floor_offset)
        mp.AddConstraint(distances[3] >= floor_offset)

        # don't leave the ground if specified
        if stay_on_ground == True:
            mp.AddConstraint(distances[0] <= np.sqrt(2)+stay_on_ground_tolerance)
            mp.AddConstraint(distances[1] <= np.sqrt(2)+stay_on_ground_tolerance)
            mp.AddConstraint(distances[2] <= np.sqrt(2)+stay_on_ground_tolerance)
            mp.AddConstraint(distances[3] <= np.sqrt(2)+stay_on_ground_tolerance)

        if fix_corner_on_ground == True:
            x_pos = get_corner_x_positions(total_x[n,:], dim)
            # make left corner on the ground
            num, loc = corner_fix
            mp.AddConstraint(distances[num] == 1.0)
            mp.AddConstraint(x_pos[num] == loc)

    # ground forces can't pull on the ground
    for n in range(N):
        force = total_f[n]
        for j in range(8):
            mp.AddConstraint(force[j] <= max_ground_force)
            mp.AddConstraint(force[j] >= 0)

    # add complimentary constraint
    for n in range(N):
        force = total_f[n]
        state = total_x[n]
        theta = state[dim]

        distances = get_corner_distances(total_x[n+1,:], dim)

        s = sin(theta)
        c = cos(theta)

        z_0 = force[0]*c + force[1]*s
        z_1 = - force[2]*s + force[3]*c
        z_2 = - force[4]*c - force[5]*s
        z_3 = force[6]*s - force[7]*c

        xy_0 = - force[0]*s + force[1]*c
        xy_1 = - force[2]*c - force[3]*s
        xy_2 = force[4]*s - force[5]*c
        xy_3 = force[6]*c + force[7]*s

        mp.AddConstraint(xy_0 <= z_0*mu)
        mp.AddConstraint(xy_0 >= -z_0*mu)
        mp.AddConstraint(xy_1 <= z_1*mu)
        mp.AddConstraint(xy_1 >= -z_1*mu)
        mp.AddConstraint(xy_2 <= z_2*mu)
        mp.AddConstraint(xy_2 >= -z_2*mu)
        mp.AddConstraint(xy_3 <= z_3*mu)
        mp.AddConstraint(xy_3 >= -z_3*mu)

        # vector_0 = force[0] * force[1]
        # vector_1 = force[2] * force[3]
        # vector_2 = force[4] * force[5]
        # vector_3 = force[6] * force[7]
        # val = np.asarray([vector_0, vector_1, vector_2, vector_3])
        # mp.AddConstraint(val.dot(distances) <= complimentarity_constraint_thresh)
        # mp.AddConstraint(val.dot(distances) >= -complimentarity_constraint_thresh)

        val_0 = np.asarray([force[0], force[2], force[4], force[6]])
        val_1 = np.asarray([force[1], force[3], force[5], force[7]])
        mp.AddConstraint(val_0.dot(distances) <= complimentarity_constraint_thresh)
        mp.AddConstraint(val_0.dot(distances) >= -complimentarity_constraint_thresh)
        mp.AddConstraint(val_1.dot(distances) <= complimentarity_constraint_thresh)
        mp.AddConstraint(val_1.dot(distances) >= -complimentarity_constraint_thresh)

    # initial state, no state error allowed
    for i in range(state_len):
        initial_state_error = x_over_time[0,i] - initial_state[i]
        mp.AddConstraint(initial_state_error == 0.0)

    # don't care about final wheel angle, so skip restriction in it
    state_indices = [i for i in range(0, state_len)]
    a = state_indices[0:state_len/2-1] + state_indices[state_len/2:]
    for i in a:
        # final
        final_state_error = x_over_time[-1,i] - final_state[i]
        mp.AddConstraint(final_state_error <= final_state_error_thresh)
        mp.AddConstraint(final_state_error >= -final_state_error_thresh)

    # add time constraint
    mp.AddConstraint(time_used[0] >= min_time)
    mp.AddConstraint(time_used[0] <= max_time)

    # add torque constraints
    for n in range(N):
        mp.AddConstraint(u_over_time[n,0] <= max_torque)
        mp.AddConstraint(u_over_time[n,0] >= -max_torque)

    # try to keep the velocity of the wheel in the correct direction
    # mp.AddLinearCost(x_over_time[:,-1].sum())

    # mp.AddLinearCost(-x_over_time[:,1].sum())
    # mp.AddLinearCost(-x_over_time[N//2,1])

    print "Number of decision vars", mp.num_vars()
    print(mp.Solve())

    trajectory = mp.GetSolution(x_over_time)
    input_trajectory = mp.GetSolution(u_over_time)
    force_trajectory = mp.GetSolution(f_over_time)
    t = mp.GetSolution(time_used)
    time_array = np.arange(0.0, t, t/(N+1))

    return trajectory, input_trajectory, force_trajectory, time_array

# compute the next state by solving implicitly for contact forces
def next_state_with_contact(state, u, dt, dim=2):
    """Return the next state after computing the force for the given timestep.
    Note that this is not guaranteed to work realtime or even every time.

    Keyword arguments:
    state -- the state of the cube
    u -- the input torque on the wheel
    dt -- the timestep size
    dim -- the dimension of the state space (default 2)
    """

    mu = 0.1 # friction force
    max_ground_force = 200
    complimentarity_constraint_thresh = 0.01

    # use SNOPT in Drake
    mp = MathematicalProgram()

    floor_offset = -0.01 # used to allow a little penitration

    x = mp.NewContinuousVariables(len(state), "x_%d" % 0)
    u_decision = mp.NewContinuousVariables(len(u), "u_%d" % 0)
    f = mp.NewContinuousVariables(8, "f_%d" % 0)

    # starting values
    for i in range(len(x)):
        mp.AddConstraint(x[i] == state[i])
    for i in range(len(u)):
        mp.AddConstraint(u_decision[i] == u[i])

    dynamic_state_next = x[:] + get_nd_dynamics(x[:], u_decision[:], f[:], dim)*dt

    # can't penitrate the floor
    distances = get_corner_distances(dynamic_state_next, dim)
    mp.AddConstraint(distances[0] >= floor_offset)
    mp.AddConstraint(distances[1] >= floor_offset)
    mp.AddConstraint(distances[2] >= floor_offset)
    mp.AddConstraint(distances[3] >= floor_offset)

    # ground forces can't pull on the ground
    for j in range(8):
        # mp.AddConstraint(f[j] <= max_ground_force)
        mp.AddConstraint(f[j] >= 0)

    # add complimentary constraint
    theta = state[dim]

    distances = get_corner_distances(state, dim)

    s = sin(theta)
    c = cos(theta)

    z_0 = f[0]*c + f[1]*s
    z_1 = - f[2]*s + f[3]*c
    z_2 = - f[4]*c - f[5]*s
    z_3 = f[6]*s - f[7]*c

    xy_0 = - f[0]*s + f[1]*c
    xy_1 = - f[2]*c - f[3]*s
    xy_2 = f[4]*s - f[5]*c
    xy_3 = f[6]*c + f[7]*s

    mp.AddConstraint(xy_0 <= z_0*mu)
    mp.AddConstraint(xy_0 >= -z_0*mu)
    mp.AddConstraint(xy_1 <= z_1*mu)
    mp.AddConstraint(xy_1 >= -z_1*mu)
    mp.AddConstraint(xy_2 <= z_2*mu)
    mp.AddConstraint(xy_2 >= -z_2*mu)
    mp.AddConstraint(xy_3 <= z_3*mu)
    mp.AddConstraint(xy_3 >= -z_3*mu)

    vector_0 = f[0] * f[1]
    vector_1 = f[2] * f[3]
    vector_2 = f[4] * f[5]
    vector_3 = f[6] * f[7]

    val = np.asarray([vector_0, vector_1, vector_2, vector_3])

    mp.AddConstraint(val.dot(distances) <= complimentarity_constraint_thresh)
    mp.AddConstraint(val.dot(distances) >= -complimentarity_constraint_thresh)

    # print "Number of decision vars", mp.num_vars()
    # print(mp.Solve())
    mp.Solve()

    f_comp = mp.GetSolution(f)
    return state + get_nd_dynamics(state, u, f_comp, dim)*dt
