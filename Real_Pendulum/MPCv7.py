from Pendulum.Model import RealInvertedPendulum, SimulatedInvertedPendulum, pend_model
from Pendulum.LQR import StabilizationLQR
from Pendulum.Controller_ACADOS import MPC
from Pendulum.Plot import Plot
from casadi import *
import matplotlib.pyplot as plt
import time


def main():

    simulation = False

    if simulation:
        x = np.array([np.pi, 0, 0, 0])
        pend = SimulatedInvertedPendulum(g=9.81, l=0.41, b=0.016, m=0.3, h=2/100, x0=x)
        stop_time = 8
        stop_iter = stop_time/pend.h
        controller_type = "MPC"
        if controller_type == 'MPC':
            goal = np.array([0, 0, 0, 0])
            contr = MPC(Tf=5, freq=250, pendulum_object=pend, xr=goal,simulation=simulation)
        else:
            contr = StabilizationLQR(pend)

    else:
        port = '/dev/ttyUSB0'                           # USB port connected to the ESP32
        baud_rate = 115200                              # Baud rate in bits per second
        stop_time = 8                                 # Stop time in seconds
        pend = RealInvertedPendulum(g=9.81, l=0.41, b=0.016, m=0.3, h=2/100, port=port, baud_rate=baud_rate)
        stop_iter = stop_time / pend.h
        controller_type = "MPC"
        if controller_type == 'MPC':
            goal = np.array([0, 0, 0, 0])
            contr = MPC(Tf=5, freq=250, pendulum_object=pend, xr=goal,port=port, baud_rate=baud_rate,simulation=simulation)
        else:
            contr = StabilizationLQR(pend)

    x_traj_meas = [pend.x.reshape(pend.nx, 1)]
    x_traj_calc = [pend.x.reshape(pend.nx, 1)]
    u_traj_calc = []
    u_traj_meas = []
    counter_iteration = 0
    new_plot = Plot(pend.h, stop_iter, simulation)
    t_loop = []
    t_ocp = []

    """Main MPC Loop"""

    while True:
        t_a = time.perf_counter()                        # Start timing block

        x = pend.measure()                               # Measuring current state

        t_ocp_a = time.perf_counter()
        u, x_k, opt_u, opt_x = contr.step(x)             # Calculating next optimal input
        t_ocp_b = time.perf_counter()

        pend.apply(u)                                    # Applying input

        x_traj_meas.append(x)                            # Measured x trajectory
        u_traj_calc.append(u)                            # u trajectory

        if simulation:
            u_traj_meas.append(u)                        # Calculated u (Since simulation doesn't involve data transfer)
        else:
            u_traj_meas.append(pend.ser.u_prev_dec)      # Applied u sent by the ESP

        counter_iteration += 1
        if stop_iter is not None:
            if counter_iteration >= stop_iter:           # Break out of loop if time exceeds stoppage time
                break

        x_traj_calc.append(x_k.reshape(4, 1))            # Calculated states
        result = contr.res                               # Controller result
        new_plot.update_plot(pend.h,
                             x_traj_meas, x_traj_calc,
                             u_traj_calc, u_traj_meas,
                             result )

        t = time.perf_counter() - t_a
        if pend.h - t >= 0:                              # Sleep for the remaining time
            time.sleep(pend.h - t)

        t_loop.append(time.perf_counter() - t_a)
        t_ocp.append(t_ocp_b - t_ocp_a)                  # Time taken to solve OCP at each iteration

    fig, ax = plt.subplots()  # plt.figure()
    tx = np.arange(0, stop_iter-1)*pend.h
    ax.plot(tx, t_ocp, label='OCP timing')
    ax.plot(tx, t_loop, label='loop timing')
    ax.legend()
    ax.set_xlabel('$t$ [s]')
    ax.set_ylabel('$t_{ocp}$ [s]')
    fig.suptitle('OCP timing')
    plt.show()


if __name__ == '__main__':
    main()
