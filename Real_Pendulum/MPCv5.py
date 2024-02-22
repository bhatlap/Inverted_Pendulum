from casadi import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import time
from statistics import mean
import scipy.linalg as lin
from scipy import signal
import serial
import struct
# import matplotlib.animation as ani
import pandas as pd
# from IPython import display


def lqr(obj, Q, R, dt):

    x = SX.sym('x',obj.nx,1)
    u = SX.sym('u',obj.nu,1)
    xdot = obj.system(x,u)
    nx = obj.nx
    nu = obj.nu

    A_fun = Function('Linearized_A',[x,u],[jacobian(xdot,x)])
    B_fun = Function('Linearized_B',[x,u],[jacobian(xdot,u)])

    xs = np.zeros([nx,1])
    us = np.zeros([nu,1])

    A = A_fun(xs,us)
    B = B_fun(xs,us)
    C = np.eye(nx)
    D = np.zeros([nx,nu])

    sys = signal.StateSpace(A,B,C,D)
    sys_d = sys.to_discrete(dt=dt)
    A_d = sys_d.A
    B_d = sys_d.B

    P = lin.solve_discrete_are(A_d,B_d,Q,R)

    return P,B_d


class PlotTrajectories:
    def __init__(self, dt, stopIter, title=" "):
        self.fig = plt.figure(num=None, figsize=(5, 10), dpi=110) # Create a figure for plotting
        self.fig.suptitle(title) # Title
        self.axes = [0, 0, 0, 0, 0] # Initialization of axis
        self.line = [0, 0, 0, 0, 0] # and line
        self.line2 = [0, 0 ,0 ,0 ,0]

        label1 = ['$\Theta$ sim', '$\dot{\Theta}$ sim', '$x$', '$\dot{x} sim$', '$u$']
        label2 = ['$\Theta$ real', '$\dot{\Theta}$ real', '$x$ real', '$\dot{x}$ real', '$u$']

        for i in range(5):
            self.axes[i] = self.fig.add_subplot(5, 1, i + 1)
            self.line[i], = self.axes[i].plot([], [], label=label1[i])
            self.line2[i], = self.axes[i].plot([], [], label=label2[i])
            self.line2[i].set_color("red")

        self.fig.tight_layout(pad=2) # Makes axis fit to data

        self.axes[0].set_ylabel('$\Theta [rad]$')
        self.axes[1].set_ylabel('$\dot{\Theta} [rad/s]$')
        self.axes[2].set_ylabel('$ p[m]$')
        self.axes[3].set_ylabel('$\dot{p} [m/s]$')
        self.axes[4].set_ylabel('$u [m/s^2]$')

        if stopIter is None:
            stopIter = 100

        # Setting x axis parameters
        for ax in self.axes:
            ax.set_xlabel("$t$ [s]")
            ax.grid()
            ax.legend()
            ax.set_xlim([0, dt*stopIter])

        # Axis limits
        self.axes[0].set_ylim([-6.5, 6.5])
        self.axes[1].set_ylim([-23, 23])
        self.axes[2].set_ylim([-0.5, 0.5])
        self.axes[3].set_ylim([-10, 10])
        self.axes[4].set_ylim([-23, 23])

        plt.show(block=False)
        self.fig.canvas.draw()

        # Lines for each state and input
        for i in range(5):
            self.axes[i].draw_artist(self.line[i])
            self.axes[i].draw_artist(self.line2[i])

        # draw the animated artist, this uses a cached renderer
        self.bg = self.fig.canvas.copy_from_bbox(self.fig.bbox)

        # show the result to the screen, this pushes the updated RGBA buffer from the
        # renderer to the GUI framework so you can see it
        self.fig.canvas.blit(self.fig.bbox)

    def update_plot(self, x_trajectory, u_trajectory, u_trajectory2,dt, extra_trajectory=" "):
        x_trajectory_conc = np.concatenate(x_trajectory, axis=1)
        u_trajectory_conc = np.array(u_trajectory)
        u_trajectory_conc2 = np.array(u_trajectory2)
        time = np.arange(0, x_trajectory_conc.shape[1]) * dt

        for i in range(4):
            self.line[i].set_data(time, x_trajectory_conc[i,:])

        self.line[4].set_data(time[:-1], u_trajectory_conc)
        self.line2[4].set_data(time[:-1], u_trajectory_conc2)

        self.fig.canvas.restore_region(self.bg)

        for i in range(5):

            self.axes[i].draw_artist(self.line[i])
            self.axes[4].draw_artist(self.line2[i])

        self.fig.canvas.blit(self.fig.bbox)
        self.fig.canvas.flush_events()


class SerialCom:
    def __init__(self, port, baud):
        self.ser = serial.Serial(port, baud)  # open serial port 115200/256000
        time.sleep(3)
        self.u_rec = 0
        self.t = time.perf_counter()
        self.u_recd = 0
        self.x = np.array([0, 0, 0, 0]).reshape(4, 1)

    def apply(self, u):

        inpuT = struct.pack('f', u)
        self.ser.write(inpuT)
        self.ser.write(b"\n")

    def measure(self):
        buffer = self.ser.read_all()
        buffer = buffer[-42:]
        print("current buffer len: ", len(buffer))
        if len(buffer) < 21:
            x = self.x
        else:
            for i in range(21, -1, -1):
                byte = buffer[i:i + 1]
                if byte == b"?":
                    print("message buffer len: ", len(buffer[i + 1:i + 21]))
                    u_recn = buffer[i + 1:i + 5]
                    if u_recn != self.u_rec:
                        print("time taken to update u ", time.perf_counter() - self.t)
                        self.t = time.perf_counter()

                    self.u_rec = u_recn
                    u_recn = struct.unpack_from('f', u_recn, offset=0)[0]
                    self.u_recd = u_recn
                    print(u_recn)

                    # Theta
                    data1 = buffer[i + 5:i + 9]
                    data1 = struct.unpack_from('f', data1, offset=0)[0]

                    # Theta_dot
                    data2 = buffer[i + 9:i + 13]
                    data2 = struct.unpack_from('f', data2, offset=0)[0]

                    # Position
                    data3 = buffer[i + 13:i + 17]
                    data3 = struct.unpack_from('f', data3, offset=0)[0]

                    # Velocity
                    data4 = buffer[i + 17:i + 21]
                    data4 = struct.unpack_from('f', data4, offset=0)[0]

                    if (abs(data1) > 1000) or (abs(data2) > 1000) or (abs(data3) > 1000) or (abs(data4) > 1000) or (abs(u_recn) > 1000):
                        print("oops")
                        print(buffer)
                        print(data1, data2, data3, data4)

                    x = np.array([data1, data2, data3, data4]).reshape(4, 1)
                    self.x = x
                    break
        return x


class OCP:
    def __init__(self, T, h, pend, lqr, xr):
        self.T = T  # Prediction Horizon in seconds
        self.h = h
        self.pend = pend  # Object of Pendulum class
        self.lqr = lqr
        self.xr = xr  # Set point

        self.N = int(self.T / self.h)
        self.Q = np.diag([1e2, 1e0, 1e2, 1e0])
        self.R = 1e1 * np.eye(pend.nu)

        self.J = 0
        self.w = 0  # Decision variables
        self.lbw = 0  # Initializing Lower bound on decision variable
        self.ubw = 0  # Initializing Upper bound on decision variable

    def build(self, verbose=True, state_constraint=True, input_constraint=True):
        """
        Function that creates and returns an optimization solver for the optimization problem defined as:

        minimize     sum_{k=0}^{N-1} (x(k)-xr).T*Q*(x(k)-xr) + u(k).T*R*u(k) + x.T*Qt*x
        x,u

        This function returns the solver object
        """
        self.J = 0
        x = SX.sym('x', self.pend.nx, 1)
        u = SX.sym('u', self.pend.nu, 1)

        X = SX.sym('X', (self.N + 1) * self.pend.nx, 1)  # Vector to store all states for N time steps
        U = SX.sym('U', self.N * self.pend.nu)  # Vector to store inputs for N time steps

        Qt, B = self.lqr(self.pend, self.Q, self.R, self.h)

        # Stage cost as a function of current state and input
        stage_cost = (x-self.xr).T @ self.Q @ (x-self.xr) + u.T @ self.R @ u
        stage_cost_fcn = Function('stage_cost', [x, u], [stage_cost])

        # Terminal cost as a function of final state
        terminal_cost = (x - self.xr).T @ Qt @ (x - self.xr)
        terminal_cost_fcn = Function('terminal_cost', [x], [terminal_cost])

        # state constraints
        lb_st = []
        ub_st = []

        # input constraints"
        lb_u = []
        ub_u = []

        G = []  # Constraint list

        # Adding initial condition to the list of constraints
        X0 = SX.sym('x_init', self.pend.nx)
        g0 = X[0:self.pend.nx] - X0
        G.append(g0)

        for k in range(self.N):

            # Extracting current, next states and current input from respective vectors

            x_k = X[k * self.pend.nx:(k + 1) * self.pend.nx, :]
            x_k_plus = X[(k + 1) * self.pend.nx:(k + 2) * self.pend.nx, :]
            u_k = U[k * self.pend.nu:(k + 1) * self.pend.nu, :]
            self.J += stage_cost_fcn(x_k, u_k)

            # Continuity constraint
            xplus = self.pend.simulate(x_k, u_k)  # Predicting the next state of the system using discretized system
            gk = xplus - x_k_plus
            G.append(gk)  # Difference between predicted state and next state in list added to constraint list

            # Updating constraints for states and inputs
            if state_constraint:
                lb_st += self.pend.x_min
                ub_st += self.pend.x_max
            else:
                lb_st += [-inf, -inf, -inf, -inf]
                ub_st += [inf, inf, inf, inf]

            if input_constraint:
                lb_u += self.pend.u_min
                ub_u += self.pend.u_max
            else:
                lb_u += [-inf]
                ub_u += [inf]

        # Terminal cost
        x_n = X[self.N * self.pend.nx:]
        self.J += terminal_cost_fcn(x_n)

        if state_constraint:
            lb_st += self.pend.x_min
            ub_st += self.pend.x_max
        else:
            lb_st += [-inf, -inf, -inf, -inf]
            ub_st += [inf, inf, inf, inf]

        # Concatenating the required vectors
        self.lbw = vertcat(*lb_st, *lb_u)
        self.ubw = vertcat(*ub_st, *ub_u)
        self.w = vertcat(X, U)

        prob = {"x": self.w, "f": self.J, "g": vertcat(*G), 'p': X0}

        #opts = {'ipopt.tol': 0.01}
        opts = {}
        if not verbose:
            opts = {**opts,
                    'ipopt.print_level': 0,
                    'print_time': 0
                    }
        else:
            opts = {**opts,
                    'ipopt.print_level': 5,
                    'print_time': 0
                    }
        opt_solver = nlpsol('solver', 'ipopt', prob, opts)

        return opt_solver

    def solve(self, start_point, init_guess, solver):
        """
        This function returns the optimal x and u trajectories after solving the optimization problem defined in the build function

        Input arguments: init_guess : Initial guess to the solver object
        start_point : State from where the optimal x trajectory will start (This is passed as a parameter to the solver to ensure continuity)

        Outputs:
            opt_x,opt_u : Optimal x and u trajectories
        """
        sol = solver(x0=init_guess, p=start_point, lbx=self.lbw, ubx=self.ubw, lbg=0, ubg=0)

        opt_x = sol['x'][:(self.N + 1) * self.pend.nx].full().reshape(self.N + 1, self.pend.nx)
        opt_u = sol['x'][(self.N + 1) * self.pend.nx:].full().reshape(self.N, self.pend.nu)

        return opt_x, opt_u


class SwingUpLQR:
    def __init__(self, lqr, pend):
        self.lqr = lqr
        self.pend = pend
        self.Q = np.diag([10, 0.1, 10, 0.1])
        self.R = 1 * np.eye(self.pend.nu)

        self.P = 0
        self.B = 0
        self.K = 0

    def step(self, x):
        self. P, self.B = self.lqr(self.pend,self.Q,self.R,self.pend.h)
        self.K = 1/self.R * self.B.T @ self.P
        u = - np.dot(self.K, x)

        if u < -23:
            u = -23
        elif u > 23:
            u = 23

        return u


class Controller:
    def __init__(self, ocp, pend, dt):

        self.ocp = ocp  # OCP object
        self.obj = pend  # System Object
        self.solver = self.ocp.build()

        self.dt = dt
        self.ocp_steps_per_mpc_iter = int(self.dt/self.obj.h)
        self.res_x_mpc = [self.obj.x]
        self.res_u_mpc = []
        self.opt_x = 0
        self.opt_u = 0
        self.ready = 0
        self.counter = 0
        self.u_k = 0

    def step(self, x):
        if self.counter % self.ocp_steps_per_mpc_iter == 0:
            self.u_k = self.calc(x)
            self.counter = 0
        u = self.u_k[self.counter]
        self.counter += 1

        return u

    def calc(self, x):

        if not self.ready:
            # First initial guess
            init_x, init_u = (np.repeat(x, self.ocp.N + 1, axis=1).T, [0] * self.ocp.N)
            init_guess = np.concatenate((init_x.flatten(), init_u[0:self.ocp.N]))

        else:
            # Warm-starting the solver
            opt_u_next = self.opt_u.flatten()[self.obj.nu:self.obj.nu * self.ocp.N]
            opt_u_next = np.concatenate((opt_u_next, self.opt_u[self.ocp.N - 1]), axis=0).reshape(self.ocp.N * self.obj.nu, 1)

            opt_x_next = self.opt_x.flatten()[self.obj.nx:self.obj.nx * (self.ocp.N + 1)].reshape(self.ocp.N * self.obj.nx, 1)
            extension = self.obj.simulate(self.opt_x[self.ocp.N], self.opt_u[self.ocp.N - 1])
            opt_x_next = np.concatenate((opt_x_next, extension), axis=0)

            init_guess = np.concatenate((opt_x_next, opt_u_next), axis=0)

        self.opt_x, self.opt_u = self.ocp.solve(x, init_guess, self.solver)

        if not self.solver.stats()["success"]:
            print('OCP not solved') # it doesn't necessarily mean that the OCP is infeasible,
                                    # just that ipopt wasn't capable of finding a solution
        else:
            print("OCP solved")

        u_k = self.opt_u[0:self.ocp_steps_per_mpc_iter]

        self.ready = 1
        return u_k


class AbstractInvertedPendulum:
    def __init__(self, g, l, b, m, h):
        self.nx = 4  # Number of states
        self.nu = 1  # Number of inputs

        self.x_min = [-np.inf, -np.inf, -0.25, -np.inf]  # State lower limits
        self.x_max = [np.inf, np.inf, 0.25, np.inf]  # State upper limits

        self.u_min = [-23]  # Input lower bound
        self.u_max = [23]  # Input upper bound

        self.h = h  # Discretisation rate

        # Constants
        self.g = g
        self.l = l
        self.m = m
        self.b = b

    def system(self, x, u):
        """
    System function that returns the next state provided the current state and input
    Provides the continuous time dynamics

    input arguments:
    a. x - Current state
    b. u - Current input
        """
        xdot = vertcat(x[1],
                       (3 * self.g / (2 * self.l)) * np.sin(x[0]) - (3 * self.b / (self.m * self.l ** 2)) * x[1] + 3 / (
                                   2 * self.l) * np.cos(x[0]) * u,
                       x[3],
                       u
                       )
        return xdot

    def simulate(self, x, u):
        """
        Function to simulate the system for one time step, given the current state and input
        """
        k1 = self.system(x, u)
        k2 = self.system(x + self.h / 2 * k1, u)
        k3 = self.system(x + self.h / 2 * k2, u)
        k4 = self.system(x + self.h * k3, u)
        xnext = x + self.h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

        return xnext


class SimulatedInvertedPendulum(AbstractInvertedPendulum):
    def __init__(self, g, l, b, m, h, x0):
        super(SimulatedInvertedPendulum, self).__init__(g, l, b, m, h)

        self.x = x0.reshape(self.nx, 1)

    def step(self, u):
        """
    Function for finding the next state using the Runge-Kutta 4 discretisation scheme
    Discrete time dynamics

    input argument:
    a. u - current input

    output:
    a. x - next state
        """
        k1 = self.system(self.x, u)
        k2 = self.system(self.x + self.h / 2 * k1, u)
        k3 = self.system(self.x + self.h / 2 * k2, u)
        k4 = self.system(self.x + self.h * k3, u)
        self.x = self.x + self.h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        self.x = self.x.full()

        return self.x

    def measure(self):
        return self.x

    def apply(self, u):

        k1 = self.system(self.x, u)
        k2 = self.system(self.x + self.h / 2 * k1, u)
        k3 = self.system(self.x + self.h / 2 * k2, u)
        k4 = self.system(self.x + self.h * k3, u)
        self.x = self.x + self.h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        self.x = self.x.full()


class RealInvertedPendulum(AbstractInvertedPendulum):
    def __init__(self, g, l, b, m, h, port, baud_rate):
        super(RealInvertedPendulum, self).__init__(g, l, b, m, h)

        self.ser = SerialCom(port, baud_rate)
        self.x = self.ser.measure()

    def step(self, u):
        self.ser.apply(u)
        time.sleep(self.h)
        self.x = self.ser.measure()

        return self.x

    def apply(self, u):
        self.ser.apply(u)

    def measure(self):
        # dig the latest state value from the buffer
        return self.ser.measure()


def main():

    simulation = True
    if simulation:
        x = np.array([np.pi, 0, 0, 0])
        pend = SimulatedInvertedPendulum(g=9.81, l=0.41, b=0.016, m=0.3, h=1/100, x0=x)
        stop_time = 8
        stop_iter = stop_time/pend.h
        if stop_iter > 600:
            stop_iter = 600

    else:
        port = 'COM4'
        baud_rate = 115200
        stop_iter = 200
        pend = RealInvertedPendulum(g=9.81, l=0.41, b=0.016, m=0.3, h=1/100, port=port, baud_rate=baud_rate)

    x_trajectory = [pend.x.reshape(pend.nx, 1)]
    u_trajectory = []
    u_trajectory2 = []
    controller_type = "mpc"

    if controller_type == "mpc":
        # Set point
        xr = np.array([0, 0, 0, 0]).reshape(pend.nx, 1)
        dt = 1/100
        ocp = OCP(1.5, pend.h, pend, lqr, xr)
        contr = Controller(ocp, pend, dt)
        start_delay = True

    else:
        contr = SwingUpLQR(lqr, pend)
        start_delay = True

    counter_iteration = 0

    new_plot = PlotTrajectories(pend.h, stop_iter)
    t_ocp = []

    x = pend.x  # Measuring initial state
    time.sleep(2)
    while True:
        # measure pendulum.meausre() function readall
        x = pend.measure()
        if start_delay:  # Delaying states by applying no input
            for k in range(100):
                x = pend.step(0)
                u = 0
            start_delay = False

        t_ocp_a = time.perf_counter()
        u = contr.step(x)  # Calculating next optimal input
        t_ocp_b = time.perf_counter()
        t_ocp.append(t_ocp_b - t_ocp_a)

        t = time.perf_counter() - t_ocp_a
        pend.apply(u)  # Applying input and measuring next state

        if pend.h - t >= 0:
            time.sleep(pend.h - t)

        print(x)
        x_trajectory.append(x)
        u_trajectory.append(u)

        if simulation:
            u_trajectory2.append(u)  # Calculated u
        else:
            u_trajectory2.append(pend.ser.u_recd)

        counter_iteration += 1
        if stop_iter is not None:
            if counter_iteration >= stop_iter:
                break

        new_plot.update_plot(x_trajectory, u_trajectory, u_trajectory2, pend.h)

    pend.apply(0)
    plt.figure()
    tx = np.arange(0, stop_iter)*pend.h
    plt.plot(tx, t_ocp)

    np.save('x_values', x_trajectory)
    np.save('u_values', u_trajectory)

    plt.show()


if __name__ == '__main__':
    main()

