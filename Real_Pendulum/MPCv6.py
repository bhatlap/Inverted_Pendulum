from casadi import *
import matplotlib.pyplot as plt
from Pendulum.Model import RealInvertedPendulum,SimulatedInvertedPendulum
import time
from Pendulum.LQR import StabilizationLQR
from multiprocessing import Process, Value, Array, Manager
from Pendulum.Controller_ACADOS import MPC
from acados_template import AcadosOcpSolver


def controller_function(pend, x, u, stop_flag, controller_type, start_delay, simulation):
    if controller_type == "mpc":
        # Set point
        controller = MPC(Tf=1, freq=100, pendulum_object=pend, simulation=simulation)
        start_delay.value = False
    else:
        controller = StabilizationLQR(pend)
        start_delay.value = False

    while True:
        t_ocp_a = time.perf_counter()
        x_reshaped = np.array([x[0], x[1], x[2], x[3]]).reshape(4, 1)
        print(x_reshaped)
        u[0] = controller.step(x_reshaped)
        print(u[0])
        t = time.perf_counter() - t_ocp_a
        if pend.h - t >= 0:
            time.sleep(pend.h - t)
        if stop_flag.value:
            break


def pendulum_function(pend, x, u, start_delay, stop_iter, stop_flag, x_trajectory, u_trajectory):
    if start_delay.value:  # Delaying states by applying no input
        for k in range(100):
            t_ocp_a = time.perf_counter()
            x[:] = pend.measure()
            pend.apply(0)  # Applying input and measuring next state
            t = time.perf_counter() - t_ocp_a
            if pend.h - t >= 0:
                time.sleep(pend.h - t)
        start_delay.value = False
    counter_iteration = 0
    while True:
        t_ocp_a = time.perf_counter()
        x[:] = pend.measure()
        pend.apply(u[0])  # Applying input and measuring next state
        t = time.perf_counter() - t_ocp_a

        x_trajectory.append(np.array([x[0], x[1], x[2], x[3]]).reshape(4, 1))
        u_trajectory.append(u[0])

        if pend.h - t >= 0:
            time.sleep(pend.h - t)

        counter_iteration += 1
        print(counter_iteration)
        if stop_iter is not None:
            if counter_iteration >= stop_iter:
                stop_flag.value=True
                break


def plotting_function(pend, x, u, stop_flag, x_trajectory, u_trajectory):
    fig, axs = plt.subplots(5, 1, sharex=True, figsize=(8, 8))
    line = [0, 0, 0, 0, 0]  # Initialization of line

    label = ['$\Theta$', '$\dot{\Theta}$', '$p$', '$\dot{p}$', '$u$']
    ylabel = ['$\Theta [rad]$', '$\dot{\Theta} [rad/s]$', '$ p[m]$', '$\dot{p} [m/s]$', '$u [m/s^2]$']

    # x_trajectory = [pend.x.reshape(pend.nx, 1)]
    # u_trajectory = []
    plt.show(block=False)

    first_time = True
    while True:
        # x_trajectory.append(np.array([x[0], x[1], x[2], x[3]]).reshape(4, 1))
        # u_trajectory.append(u[0])
        x_list = list(x_trajectory)
        u_list = list(u_trajectory)

        # x_conc=np.array(x_list)
        x_trajectory_conc = np.concatenate(x_list, axis=1)
        u_trajectory_conc = np.array(u_list)
        #
        # print(x_conc)
        # print(u_list)

        time = np.arange(0, x_trajectory_conc.shape[1]) * pend.h
        print(time[-1])

        for i in range(4):
            line[i], = axs[i].plot(time, x_trajectory_conc[i, :], label=label[i])
            line[i].set_color("blue")
        line[4], = axs[4].plot(time[:-1], u_trajectory_conc[0:x_trajectory_conc.shape[1]-1], label=label[4])
        line[4].set_color("blue")

        if first_time:
            for i in range(5):
                axs[i].set_xlabel("$t$ [s]")
                axs[i].set_ylabel(ylabel[i])
                axs[i].grid()
                axs[i].legend()
            first_time = False

        fig.tight_layout(pad=2)
        plt.pause(0.000001)
        if stop_flag.value:
            plt.show()
            break


def main():

    simulation = False
    if simulation:
        x = np.array([np.pi, 0, 0, 0])
        pend = SimulatedInvertedPendulum(g=9.81, l=0.41, b=0.016, m=0.3, h=0.01, x0=x)
        stop_time = 16
        stop_iter = stop_time/pend.h
        # if stop_iter > 200:
        #     stop_iter = 200

    else:
        port = '/dev/ttyUSB0'
        baud_rate = 115200
        stop_iter = None
        pend = RealInvertedPendulum(g=9.81, l=0.41, b=0.016, m=0.3, h=1/100, port=port, baud_rate=baud_rate)

    controller_type = "mpc"

    u_shared = Array('f', range(1))
    x_shared = Array('f', range(4))
    stop_flag = Value('b', False)
    start_delay = Value('b', False)
    x_shared[:] = pend.x  # Measuring initial state

    manager = Manager()
    x_lst = manager.list()
    u_lst = manager.list()
    x_lst.append(pend.x.reshape(pend.nx, 1))

    time.sleep(2)
    p1 = Process(target=controller_function, args=(pend, x_shared, u_shared, stop_flag, controller_type, start_delay, simulation))
    p2 = Process(target=pendulum_function, args=(pend, x_shared, u_shared, start_delay, stop_iter, stop_flag, x_lst, u_lst))
    p3 = Process(target=plotting_function, args=(pend, x_shared, u_shared, stop_flag, x_lst, u_lst))
    p1.start()
    p2.start()
    p3.start()

    p1.join()
    p2.join()
    p3.join()
    print(1)

if  __name__== '__main__':
    main()

