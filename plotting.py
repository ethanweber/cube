import matplotlib.pyplot as plt

def xy_plot(traj, time_array, dim=2):
    """Creates a plot of the x,y position over time.

    Keyword arguments:
    traj -- trajectory array
    time_array -- time corresponding to state
    dim -- state dimension (default 2)
    """
    x_traj = traj[:,0]
    y_traj = traj[:,1]

    fig, ax = plt.subplots()

    ax.plot(time_array[:len(traj)], x_traj, 'b--', label='x')
    ax.plot(time_array[:len(traj)], y_traj, 'r', label='y')

    legend = ax.legend(shadow=True, fontsize='x-large')
    legend.get_frame().set_facecolor('#FFFFFF')

    plt.title("x and y vs. time")
    plt.xlabel("time (sec)")
    plt.show()

def x_plot(traj, time_array, dim=2):
    """Creates a plot of the x position over time.

    Keyword arguments:
    traj -- trajectory array
    time_array -- time corresponding to state
    dim -- state dimension (default 2)
    """
    x_traj = traj[:,0]

    fig, ax = plt.subplots()

    ax.plot(time_array[:len(traj)], x_traj, 'b--', label='x')

    legend = ax.legend(shadow=True, fontsize='x-large')
    legend.get_frame().set_facecolor('#FFFFFF')

    plt.title("x vs. time")
    plt.xlabel("time (sec)")
    plt.show()

def y_plot(traj, time_array, dim=2):
    """Creates a plot of the y position over time.

    Keyword arguments:
    traj -- trajectory array
    time_array -- time corresponding to state
    dim -- state dimension (default 2)
    """
    y_traj = traj[:,1]

    fig, ax = plt.subplots()

    ax.plot(time_array[:len(traj)], y_traj, 'r', label='y')

    legend = ax.legend(shadow=True, fontsize='x-large')
    legend.get_frame().set_facecolor('#FFFFFF')

    plt.title("y vs. time")
    plt.xlabel("time (sec)")
    plt.show()

def theta_plot(traj, time_array, dim=2):
    """Creates a theta plot.

    Keyword arguments:
    traj -- trajectory array
    time_array -- time corresponding to state
    dim -- state dimension (default 2)
    """
    theta_traj = traj[:,dim]

    fig, ax = plt.subplots()

    ax.plot(time_array[:len(traj)], theta_traj, 'b--', label='theta')

    legend = ax.legend(shadow=True, fontsize='x-large')
    legend.get_frame().set_facecolor('#FFFFFF')

    plt.title("theta vs. time")
    plt.xlabel("time (sec)")
    plt.show()

def alpha_plot(traj, time_array, dim=2):
    """Creates a alpha plot.

    Keyword arguments:
    traj -- trajectory array
    time_array -- time corresponding to state
    dim -- state dimension (default 2)
    """
    alpha_index = len(traj[0])/2 - 1
    alpha_traj = traj[:,alpha_index]

    fig, ax = plt.subplots()

    ax.plot(time_array[:len(traj)], alpha_traj, 'b--', label='alpha')

    legend = ax.legend(shadow=True, fontsize='x-large')
    legend.get_frame().set_facecolor('#FFFFFF')

    plt.title("alpha vs. time")
    plt.xlabel("time (sec)")
    plt.show()


def plot_states(traj, time_array, dim=2):

    x_traj = traj[:,0]
    y_traj = traj[:,1]

    theta_traj = traj[:,dim]

    alpha_index = len(traj[0])/2 - 1
    alpha_traj = traj[:,alpha_index]

    f, axarr = plt.subplots(2, 2, figsize=(10,5))

    axarr[0, 0].plot(time_array[:len(traj)], x_traj, 'b--', label='x')
    axarr[0, 0].legend(shadow=True, fontsize='x-large').get_frame().set_facecolor('#FFFFFF')

    axarr[0, 1].plot(time_array[:len(traj)], y_traj, 'b--', label='y')
    axarr[0, 1].legend(shadow=True, fontsize='x-large').get_frame().set_facecolor('#FFFFFF')

    axarr[1, 0].plot(time_array[:len(traj)], theta_traj, 'b--', label='theta')
    axarr[1, 0].legend(shadow=True, fontsize='x-large').get_frame().set_facecolor('#FFFFFF')

    axarr[1, 1].plot(time_array[:len(traj)], alpha_traj, 'b--', label='alpha')
    axarr[1, 1].legend(shadow=True, fontsize='x-large').get_frame().set_facecolor('#FFFFFF')

    plt.suptitle('states vs. time (sec)')
    plt.show()

def plot_vels(traj, time_array, dim=2):

    hl = len(traj[0]) / 2

    x_dot_traj = traj[:,0+hl]
    y_dot_traj = traj[:,1+hl]
    theta_dot_traj = traj[:,dim+hl]
    alpha_dot_traj = traj[:,-1]

    f, axarr = plt.subplots(2, 2, figsize=(10,5))

    axarr[0, 0].plot(time_array[:len(traj)], x_dot_traj, 'b--', label='x_dot')
    axarr[0, 0].legend(shadow=True, fontsize='x-large').get_frame().set_facecolor('#FFFFFF')

    axarr[0, 1].plot(time_array[:len(traj)], y_dot_traj, 'b--', label='y_dot')
    axarr[0, 1].legend(shadow=True, fontsize='x-large').get_frame().set_facecolor('#FFFFFF')

    axarr[1, 0].plot(time_array[:len(traj)], theta_dot_traj, 'b--', label='theta_dot')
    axarr[1, 0].legend(shadow=True, fontsize='x-large').get_frame().set_facecolor('#FFFFFF')

    axarr[1, 1].plot(time_array[:len(traj)], alpha_dot_traj, 'b--', label='alpha_dot')
    axarr[1, 1].legend(shadow=True, fontsize='x-large').get_frame().set_facecolor('#FFFFFF')

    plt.suptitle('velocities vs. time (sec)')
    plt.show()

def input_plot(input_traj, time_array, dim=2):
    """Creates a plot of input torque on the wheel.

    Keyword arguments:
    input_traj -- trajectory array
    time_array -- time corresponding to state
    dim -- state dimension (default 2)
    """

    fig, ax = plt.subplots()

    ax.plot(time_array[:len(input_traj)], input_traj[:], 'b--', label='torque')

    legend = ax.legend(shadow=True, fontsize='x-large')
    legend.get_frame().set_facecolor('#FFFFFF')

    plt.title("torque vs. time (sec)")
    plt.show()

def ground_force_plot(force_traj, time_array, dim=2, side="left"):

    fig, ax = plt.subplots()

    ax.plot(time_array[:len(force_traj)], force_traj[:,0], 'c--', label='0')
    ax.plot(time_array[:len(force_traj)], force_traj[:,1], 'c', label='1')

    ax.plot(time_array[:len(force_traj)], force_traj[:,2], 'g--', label='2')
    ax.plot(time_array[:len(force_traj)], force_traj[:,3], 'g', label='3')

    ax.plot(time_array[:len(force_traj)], force_traj[:,4], 'b--', label='4')
    ax.plot(time_array[:len(force_traj)], force_traj[:,5], 'b', label='5')

    ax.plot(time_array[:len(force_traj)], force_traj[:,6], 'o--', label='6')
    ax.plot(time_array[:len(force_traj)], force_traj[:,7], 'o', label='7')

    legend = ax.legend(shadow=False, fontsize='large', loc='upper {}'.format(side))
    legend.get_frame().set_facecolor('#FFFFFF')

    plt.title("contact force vs. time (sec)")
    plt.show()


    # plt.savefig("swing_up_states.png")
