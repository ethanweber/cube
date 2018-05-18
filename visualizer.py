import meshcat
import meshcat.geometry as g
import meshcat.transformations as tf

class MeshcatCubeVisualizer:
    def __init__(self):
        self.vis = meshcat.Visualizer()

        self.mblock = None # to be drawn later if wanted

        self.cube = self.vis["cube"]
        self.pivot = self.cube["pivot"]
        self.wheel = self.pivot["wheel"]

        # create and draw the cube
        self.cube_dim = [1.0,1.0,1.0] # x,y,z
        self.cube.set_object(g.Box(self.cube_dim))

        # pivot and wheel
        self.pivot.set_transform(tf.translation_matrix([0,0,0])) # set location of pole
        wheel_dim = [1.5,.5,.5] # x,y,z
        self.wheel.set_object(g.Box(wheel_dim))

        self.initialize()

    def draw_transformation(self, state, dim=2.0):
        nd = len(state)
        state = list(state)
        origin = state[0:3] # so you can edit
        origin[0] = 0.0
        origin[1] = state[0]
        origin[2] = state[1] + self.cube_dim[2]/2.0
        theta = state[dim]
        wheel_angle = state[len(state)/2 - 1]
        temp = tf.rotation_matrix(theta,[1,0,0]) # assume rotate about y
        temp[0:3, -1] = tf.translation_from_matrix(tf.translation_matrix(origin))
        self.cube.set_transform(temp)
        self.wheel.set_transform(tf.rotation_matrix(wheel_angle,[1,0,0])) # rotate the pole

    def initialize(self):
        # set the initial state in 2d
        x = 0.0
        y = 0.0
        x_dot = 0.0
        y_dot = 0.0
        theta = 0.0
        theta_dot = 0.0
        # state of the flywheel
        alpha = 0.0
        alpha_dot = 0.0

        state_initial = (x,y,theta,alpha,x_dot,y_dot,theta_dot,alpha_dot)
        self.draw_transformation(state_initial, 2)

    def draw_mblock(self):
        self.mblock = self.vis["mblock"]
        # create and draw the mblock
        dim = [3,1,1]
        self.mblock.set_object(g.Box(dim))
        temp = tf.translation_matrix([0,1,0.5])
        self.mblock.set_transform(temp)
