# define thresholds throughout
import numpy as np
from pydrake.math import sin, cos

def get_corner_distances(state, dim=2, linearize_theta=None):
    """Return distances to ground in the corner order [0,1,2,3] as a numpy array.

    Keyword arguments:
    state -- the state of the cube
    dim -- the dimension of the state space (default 2)
    """

    y = state[1]
    if linearize_theta is not None:
        theta = linearize_theta
    else:
        theta = state[dim]

    offset = .5*np.sqrt(2)*sin(np.pi/4.0+theta)
    val = sin(theta)

    dist_0 = y - offset
    dist_1 = dist_0 + val
    dist_2 = y + offset
    dist_3 = dist_2 - val

    # add .5 for the offset
    dist_0 += .5
    dist_1 += .5
    dist_2 += .5
    dist_3 += .5

    return np.asarray([dist_0, dist_1, dist_2, dist_3])

def get_corner_x_positions(state, dim=2, linearize_theta=None):
    """Return x position of the corners in the order [0,1,2,3] as a numpy array.

    Keyword arguments:
    state -- the state of the cube
    dim -- the dimension of the state space (default 2)
    """

    x = state[0]
    if linearize_theta is not None:
        theta = linearize_theta
    else:
        theta = state[dim]

    offset = .5*np.sqrt(2)*cos(np.pi/4.0+theta)
    val = cos(theta)

    pos_0 = x - offset
    pos_1 = pos_0 + val
    pos_2 = x + offset
    pos_3 = pos_2 - val

    return np.asarray([pos_0, pos_1, pos_2, pos_3])

floor_offset = -.01 # used to allow a little penitration
def add_floor_constraint(mp, state, dim=2):
    distances = get_corner_distances(state[:], dim)
    mp.AddConstraint(distances[0] >= floor_offset)
    mp.AddConstraint(distances[1] >= floor_offset)
    mp.AddConstraint(distances[2] >= floor_offset)
    mp.AddConstraint(distances[3] >= floor_offset)


stay_on_ground_tolerance = 0.1 # tolerance for leaving the ground if contstraint is used
def stay_on_ground(mp, state, dim=2):
    # don't leave the ground if specified
    distances = get_corner_distances(state[:], dim)
    mp.AddConstraint(distances[0] <= np.sqrt(2)+stay_on_ground_tolerance)
    mp.AddConstraint(distances[1] <= np.sqrt(2)+stay_on_ground_tolerance)
    mp.AddConstraint(distances[2] <= np.sqrt(2)+stay_on_ground_tolerance)
    mp.AddConstraint(distances[3] <= np.sqrt(2)+stay_on_ground_tolerance)


def fix_corner_to_ground(mp, state, corner_index=0, x_coord=-0.5, dim=2):
    distances = get_corner_distances(state[:], dim)
    # make left corner on the ground in specified position
    x_pos = get_corner_x_positions(state, dim)
    mp.AddConstraint(distances[corner_index] == 0.0)
    mp.AddConstraint(x_pos[corner_index] == x_coord)

def add_corner_cost(mp, state, corner_index=0, x_coord=-0.5, dim=2, linearize_theta=None):
    # make left corner on the ground in specified position with quadratic cost

    distances = get_corner_distances(state[:], dim, linearize_theta)
    mp.AddQuadraticCost(distances[corner_index]*distances[corner_index])

    x_pos = get_corner_x_positions(state, dim, linearize_theta)
    x_diff = x_pos[corner_index]  - x_coord
    # mp.AddQuadraticCost(x_diff*x_diff)


max_ground_force = 1000
def dont_pull_on_ground(mp, force, dim=2):
    for j in range(len(force)):
        mp.AddConstraint(force[j] <= max_ground_force)
        mp.AddConstraint(force[j] >= 0)


complimentarity_constraint_thresh = 0.1
mu = 10.0 # friction force
def complimentarity_constraint(mp, state, force, dim=2):
    theta = state[dim]

    distances = get_corner_distances(state[:], dim)

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


# set the entire start state
def set_initial_state(mp, x, current_state, dim):
    for i in range(len(x)):
        mp.AddConstraint(x[i] == current_state[i])

def bound_abs_value(mp, var, thresh):
    for i in range(len(var)):
        mp.AddConstraint(var[i] <= thresh)
        mp.AddConstraint(var[i] >= -thresh)
