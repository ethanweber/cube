import numpy as np
from pydrake.math import sin, cos
from pydrake.all import (SignalLogger, CompliantMaterial, ConstantVectorSource, DirectCollocation, DiagramBuilder, FloatingBaseType,
                         PiecewisePolynomial, RigidBodyTree, RigidBodyPlant,
                         SolutionResult, AddModelInstancesFromSdfString,
                         MathematicalProgram, Simulator, BasicVector, AddFlatTerrainToWorld)

def get_nd_state(state, dim=2):
    """Return higher dimensional state given a state in 2D.

    Keyword arguments:
    state -- the state of the cube in 2D
    dim -- the desired dimension (default 2)
    """
    assert len(state) == 8 # make sure input is 2D
    x,y,theta,alpha,x_dot,y_dot,theta_dot,alpha_dot = state
    if dim > 2:
        nd_state = [x] + [y] + [0]*(dim-2)
        nd_state += [theta] + [0]*(dim-1)
        nd_state += [alpha]
        nd_state += [x_dot] + [y_dot] + [0]*(dim-2)
        nd_state += [theta_dot] + [0]*(dim-1)
        nd_state += [alpha_dot]
        return nd_state
    return state

# mass
m_c = 0.5 # cube
m_w = 0.5 # wheel
# friction
F_c = 0.5
F_w = 0.5
# moments of inertia
I_c = 1.0
I_w = 1.0

def get_nd_dynamics(state, u, force, dim=2, linearize_theta=None):
    """Return state space dyanmics in n dimensions.

    Keyword arguments:
    state -- the state of the cube
    u -- the input torque on the wheel
    force -- the 8 forces acting on the wheel corners
    dim -- the dimension of the state space (default 2)
    linearize_theta -- this is used only because the qp solver needed it this way
    """
    # half state length
    hl = len(state) / 2

    m_t = m_c + m_w # total mass
    I_t = I_c + I_w # total inertia

    # gravity
    g = 9.81

    # unpack the states
    # x = state[0]
    # y = state[1]
    if linearize_theta is not None:
        theta = linearize_theta
    else:
        theta = state[dim]
    # alpha = state[hl-1]
    # xdot = state[0+hl]
    # ydot = state[1+hl]
    theta_dot = state[dim+hl]
    alpha_dot = state[-1]

    # derivative vector
    derivs = np.zeros_like(state)
    derivs[0:hl] = state[hl:] # set velocities

    # ballistic dynamics
    derivs[0+hl] = (force[1] - force[2] + force[6] - force[5])*cos(theta) - (force[0] + force[3] - force[4] - force[7])*sin(theta) # forces along x
    derivs[1+hl] = (force[1] - force[2] + force[6] - force[5])*sin(theta) + (force[0] + force[3] - force[4] - force[7])*cos(theta) - g*m_t  # forces in y direction

    # cube angle acceleration
    derivs[dim + hl] = (-u[0] + F_w*alpha_dot - F_c*theta_dot)/I_c + (-force[0]+force[1]-force[2]+force[3]-force[4]+force[5]-force[6]+force[7])*.5

    # wheel acceleration
    derivs[-1] = (u[0]*I_t + F_c*theta_dot*I_w - F_w*alpha_dot*I_t)/(I_w*I_c)

    return derivs
