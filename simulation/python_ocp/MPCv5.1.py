import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from casadi import *
import time
from statistics import mean
import scipy.linalg as lin
from scipy import signal
import serial
import struct
import matplotlib.animation as ani


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

    return P


def plot_trajectories(x, u, dt, title=""):
    time = np.arange(0, x.shape[1]) * dt
    fig, axs = plt.subplots(5, 1, sharex=True)
    fig.suptitle(title)
    axs[0].plot(time, x[0,:],label='$\Theta$')
    axs[1].plot(time, x[1,:],label='$\dot{\Theta}$')
    axs[2].plot(time, x[2,:],label='$x$')
    axs[3].plot(time, x[3,:],label='$\dot{x}$')
    axs[4].plot(time[:-1],u,label='$u$') #

    for ax in axs:
        ax.set_xlabel("$t$ [s]")
        ax.grid()
        ax.legend()
    print(f"theta at end of {title} = {x[0,-1]} ")


class SerialCom:
    def __init__(self, port, baud):
        self.ser = serial.Serial(port, baud)  # open serial port 115200/256000
        self.u_rec = 0
        self.t = time.perf_counter()

    def apply(self, u):
        # print('')
        inpuT = struct.pack('f', u)
        # print('input:')
        # print(inpuT.hex())
        self.ser.write(inpuT)
        self.ser.write(b"\n")

    def measure(self):
        byt = self.ser.read()
        while byt!=b"?":
            byt=self.ser.read()
            # print(byt)
        else:
            u_recn = self.ser.read(4)
            if u_recn != self.u_rec:
                print("time toö took to update u ", time.perf_counter()-self.t)
                self.t=time.perf_counter()
            # print(data1.hex())
            self.u_rec = u_recn
            u_recn = struct.unpack_from('f', u_recn, offset=0)[0]
            # print('u receive:')
            print(u_recn)

            data1 = self.ser.read(4)
            # print(data1.hex())
            data1 = struct.unpack_from('f', data1, offset=0)[0]
            # print(data1)

            data2 = self.ser.read(4)
            # print(data2.hex())
            data2 = struct.unpack_from('f', data2, offset=0)[0]

            # print('x:')
            # print(data1)
            data3 = self.ser.read(4)
            # print(data3.hex())
            data3 = struct.unpack_from('f', data3, offset=0)[0]

            data4 = self.ser.read(4)
            # print(data4.hex())
            data4 = struct.unpack_from('f', data4, offset=0)[0]
            x = np.array([data1, data2, data3, data4]).reshape(4,1)

            # elapsed = time.time() - self.t
            # time.sleep(0.5)
            # time.sleep(0.08 - elapsed)
            # print('freq:')
            # print(1 / (time.time() - self.t))
        return x


class OCP:
    def __init__(self, T, h, pend, lqr):
        self.T = T  # Prediction Horizon in seconds
        self.h = h
        self.pend = pend  # Object of Pendulum class
        self.lqr = lqr

        self.N = int(self.T / self.h)
        self.Q = np.diag([1e2, 1, 100, 1])
        self.R = 30 * np.eye(pend.nu)
        self.J = 0
        self.w = 0  # Decision variables
        self.lbw = 0  # Initializing Lower bound on decision variable
        self.ubw = 0  # Initializing Upper bound on decision variable

    def build(self, verbose=False, state_constraint=True, input_constraint=True):
        """
        Function that creates and returns an optimization solver for the optimization problem defined as:

        minimize     sum_{k=0}^{N-1} x(k).T*Q*x(k) + u(k).T*R*u(k) + x.T*Qt*x
        x,u

        This function returns the solver object
        """
        self.J = 0
        x = SX.sym('x', self.pend.nx, 1)
        u = SX.sym('u', self.pend.nu, 1)

        X = SX.sym('X', (self.N + 1) * self.pend.nx, 1)  # Vector to store all states for N time steps
        U = SX.sym('U', self.N * self.pend.nu)  # Vector to store inputs for N time steps

        # Qt = 100*self.Q

        Qt = self.lqr(self.pend, self.Q, self.R, self.h)

        # Stage cost as a function of current state and input
        stage_cost = x.T @ self.Q @ x + u.T @ self.R @ u
        stage_cost_fcn = Function('stage_cost', [x, u], [stage_cost])

        # Terminal cost as a function of final state
        terminal_cost = x.T @ Qt @ x
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

            # Stage cost
            self.J += stage_cost_fcn(x_k, u_k)

            # Continuity constraint
            xplus = self.pend.simulate(x_k, u_k)  # Predicting the next state of the system using discretized system
            gk = xplus - x_k_plus
            G.append(gk)  # Difference between predicted state and next state in list added to constraint list

            # Updating constraints for states and inputs
            if state_constraint:
                lb_st += [-inf, -inf, -0.15, -inf]
                ub_st += [inf, inf, 0.15, inf]
            else:
                lb_st += [-inf, -inf, -inf, -inf]
                ub_st += [inf, inf, inf, inf]

            if input_constraint:
                lb_u += [-15]
                ub_u += [15]
            else:
                lb_u += [-inf]
                ub_u += [inf]

        # Terminal cost
        x_n = X[self.N * self.pend.nx:]
        self.J += terminal_cost_fcn(x_n)

        if state_constraint:
            lb_st += [-inf, -inf, -0.15, -inf]
            ub_st += [inf, inf, 0.15, inf]
        else:
            lb_st += [-inf, -inf, -inf, -inf]
            ub_st += [inf, inf, inf, inf]

        # Concatenating the required vectors
        self.lbw = vertcat(*lb_st, *lb_u)
        self.ubw = vertcat(*ub_st, *ub_u)
        self.w = vertcat(X, U)

        prob = {"x": self.w, "f": self.J, "g": vertcat(*G), 'p': X0}

        opts = {'ipopt.tol': 1e-2}

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
    def __init__(self):
        self.K = np.array([32.7807, 4.7813, -2.9829, -3.7504])

    def step(self, x):
        return - np.dot(self.K, x)


class Controller:
    def __init__(self, ocp, pend, dt): # Include serial object to simulate arduino comms
        self.ocp = ocp  # OCP object
        self.obj = pend  # System Object
        self.solver = self.ocp.build()


        #self.ser = ser # Serial communication object

        self.dt = dt
        self.ocp_steps_per_mpc_iter = int(self.dt/self.obj.h)
        self.res_x_mpc = [self.obj.x]
        self.res_u_mpc = []
        self.opt_x = 0
        self.opt_u = 0
        self.ready_ = 0
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

        if not self.ready_:
            # Provide first initial guess
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

        self.opt_x,self.opt_u = self.ocp.solve(x,init_guess,self.solver)
        u_k = self.opt_u[0:self.ocp_steps_per_mpc_iter]

        self.ready_ = 1
        return u_k


class AbstractInvertedPendulum:
    def __init__(self, g, l, b, m, h):
        self.nx = 4  # Number of states
        self.nu = 1  # Number of inputs

        self.x_min = [-np.inf, -np.inf, -0.3, -np.inf]  # State lower limits
        self.x_max = [np.inf, np.inf, 0.3, np.inf]  # State upper limits

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


class RealInvertedPendulum(AbstractInvertedPendulum):
    def __init__(self, g, l, b, m, h, port, baud_rate):
        super(RealInvertedPendulum, self).__init__(g, l, b, m, h)

        self.ser = SerialCom(port, baud_rate)
        self.x = self.ser.measure()

    def step(self, u):
        self.ser.apply(u)
        self.x = self.ser.measure()

        return self.x

    def measure(self):
        # dig the latest state value from the buffer
        return self.ser.measure()


def main():

    simulation = False
    if simulation:
        x = np.array([np.pi, 0, 0, 0])
        pend = SimulatedInvertedPendulum(g=9.81, l=0.41, b=0.016, m=0.3, h=1/12.5, x0=x)
        stop_iter = 100
    else:
        port = 'COM4'
        baud_rate = 115200
        stop_iter = None
        pend = RealInvertedPendulum(g=9.81, l=0.41, b=0.016, m=0.3, h=1/12.5, port=port, baud_rate=baud_rate)

    x_trajectory = [pend.x.reshape(pend.nx,1)]
    u_trajectory = []

    controller_type = "mpc"

    if controller_type == "mpc":
        ocp = OCP(2, pend.h, pend, lqr)
        contr = Controller(ocp, pend, 1/12.5)
    else:
        contr = SwingUpLQR()
    counter_iteration = 0

    # stop_iter =50
    while True:
        #x = pend.measure() # getting the latest states from the receive value
        x = pend.x
        u = contr.step(x)

        x = pend.step(u)

        x_trajectory.append(x)
        u_trajectory.append(u)

        counter_iteration += 1
        if stop_iter is not None:
            if counter_iteration >= stop_iter:
                break

        print(counter_iteration)
        print(x)
        print(u)

        # sleep
        # time.sleep(1)

    x_trajectory = np.concatenate(x_trajectory, axis=1)
    u_trajectory = np.array(u_trajectory)

    plot_trajectories(x_trajectory, u_trajectory, pend.h, title="MPC")
    plt.show()


if __name__ == '__main__':
    main()