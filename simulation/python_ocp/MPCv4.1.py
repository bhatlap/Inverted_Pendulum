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
    fig, axs = plt.subplots(5, 1)
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
        self.t = 0

    def send(self, u):
        self.t = time.time()
        print('')
        inpuT = struct.pack('f', u)
        # print('input:')
        # print(inpuT.hex())
        self.ser.write(inpuT)
        self.ser.write(b"\n")

    def receive(self):
        data1 = self.ser.read(4)
        # print(data1.hex())
        data1 = struct.unpack_from('f', data1, offset=0)[0]
        # print('u receive:')
        # print(data1)

        data1 = self.ser.read(4)
        print(data1.hex())
        data1 = struct.unpack_from('f', data1, offset=0)[0]
        # print(data1)

        data2 = self.ser.read(4)
        print(data2.hex())
        data2 = struct.unpack_from('f', data2, offset=0)[0]

        # print('x:')
        # print(data1)
        data3 = self.ser.read(4)
        print(data3.hex())
        data3 = struct.unpack_from('f', data3, offset=0)[0]

        data4 = self.ser.read(4)
        print(data4.hex())
        data4 = struct.unpack_from('f', data4, offset=0)[0]
        x = np.array([data1, data2, data3, data4])

        # elapsed = time.time() - self.t
        # time.sleep(0.5)
        # time.sleep(0.08 - elapsed)
        # print('freq:')
        # print(1 / (time.time() - self.t))
        return x


class OCP:
    def __init__(self,T,h,pend,lqr):
        self.T = T # Prediction Horizon in seconds
        self.h = h
        self.pend = pend # Object of Pendulum class
        self.lqr = lqr

        self.N = int(self.T/self.h)
        self.Q = np.diag([1e2, 1, 100, 1])
        self.R = 1*np.eye(pend.nu)
        self.J = 0
        self.w = 0 # Decision variables
        self.lbw = 0 # Initializing Lower bound on decision variable
        self.ubw = 0 # Initializing Upper bound on decision variable

    def build(self, verbose=False, state_constraint=True, input_constraint=True):
        """
        Function that creates and returns an optimization solver for the optimization problem defined as:
    
        minimize     sum_{k=0}^{N-1} x(k).T*Q*x(k) + u(k).T*R*u(k) + x.T*Qt*x
        x,u

        This function returns the solver object 
        """
        self.J = 0
        x = SX.sym('x', self.pend.nx,1)
        u = SX.sym('u', self.pend.nu,1)

        X = SX.sym('X',(self.N+1)*self.pend.nx,1) # Vector to store all states for N time steps
        U = SX.sym('U',self.N*self.pend.nu) # Vector to store inputs for N time steps

        #Qt = 100*self.Q

        Qt = self.lqr(self.pend,self.Q,self.R,self.h)

        # Stage cost as a function of current state and input
        stage_cost = x.T @ self.Q @ x+ u.T @ self.R @ u
        stage_cost_fcn = Function('stage_cost',[x,u],[stage_cost])

        # Terminal cost as a function of final state
        terminal_cost = x.T @ Qt @ x
        terminal_cost_fcn = Function('terminal_cost',[x],[terminal_cost])

        # state constraints
        lb_st = []
        ub_st = []
    
        # input constraints"
        lb_u = []
        ub_u = []

        G =[] # Constraint list

        # Adding initial condition to the list of constraints
        X0 = SX.sym('x_init', self.pend.nx)
        g0 =X[0:self.pend.nx]-X0
        G.append(g0)

        for k in range(self.N):

            # Extracting current, next states and current input from respective vectors

            x_k = X[k*self.pend.nx:(k+1)*self.pend.nx,:]
            x_k_plus = X[(k+1)*self.pend.nx:(k+2)*self.pend.nx,:]
            u_k = U[k*self.pend.nu:(k+1)*self.pend.nu,:] 
        
            # Stage cost
            self.J += stage_cost_fcn(x_k, u_k)
        
            # Continuity constraint
            xplus = self.pend.simulate(x_k,u_k)  # Predicting the next state of the system using discretized system
            gk = xplus-x_k_plus
            G.append(gk) # Difference between predicted state and next state in list added to constraint list
        
            # Updating constraints for states and inputs
            if state_constraint:
                lb_st += [-inf, -inf, -0.3, -inf]
                ub_st += [inf, inf, 0.3, inf]
            else:
                lb_st += [-inf, -inf, -inf, -inf]
                ub_st += [inf, inf, inf, inf]

            if input_constraint:
                lb_u += [-23]
                ub_u += [23]
            else:
                lb_u += [-inf]
                ub_u += [inf]
            
        # Terminal cost
        x_n = X[self.N*self.pend.nx:]
        self.J += terminal_cost_fcn(x_n)
    
        if state_constraint:
            lb_st += [-inf, -inf, -0.3, -inf]
            ub_st += [inf, inf, 0.3, inf]
        else:
            lb_st += [-inf, -inf, -inf, -inf]
            ub_st += [inf, inf, inf, inf]

        # Concatenating the required vectors
        self.lbw = vertcat(*lb_st, *lb_u)
        self.ubw = vertcat(*ub_st, *ub_u)
        self.w = vertcat(X,U)

        prob = {"x": self.w,"f": self.J,"g": vertcat(*G), 'p': X0}

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
        opt_solver = nlpsol('solver','ipopt',prob, opts)

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

        opt_x = sol['x'][:(self.N+1)*self.pend.nx].full().reshape(self.N+1, self.pend.nx)
        opt_u = sol['x'][(self.N+1)*self.pend.nx:].full().reshape(self.N, self.pend.nu)

        return opt_x, opt_u


class Controller:
    def __init__(self, ocp, pend, dt): # Include serial object to simulate arduino comms
        self.ocp = ocp  # OCP object
        self.obj = pend  # System Object

        #self.ser = ser # Serial communication object
        self.solver = self.ocp.build()

        self.dt = dt
        self.res_x_mpc = [self.obj.x]
        self.res_u_mpc = []

    def mpc_sim(self, T_sim):

        N_sim = int(T_sim/self.dt)
        ocp_steps_per_mpc_iter = int(self.dt/self.obj.h)

        for i in range(N_sim):

            if i == 0:
                # Provide first initial guess
                init_x, init_u = (np.repeat(self.obj.x, self.ocp.N + 1, axis=1).T, [0] * self.ocp.N)
                init_guess = np.concatenate((init_x.flatten(), init_u[0:self.ocp.N]))

            elif i > 0:
                # Warmstart the solver
                opt_u_next = opt_u.flatten()[self.obj.nu:self.obj.nu * self.ocp.N]
                opt_u_next = np.concatenate((opt_u_next, opt_u[self.ocp.N - 1]), axis=0).reshape(self.ocp.N*self.obj.nu,1)

                opt_x_next = opt_x.flatten()[self.obj.nx:self.obj.nx *(self.ocp.N + 1)].reshape(self.ocp.N*self.obj.nx, 1)
                extension = self.obj.simulate(opt_x[self.ocp.N], opt_u[self.ocp.N - 1])
                opt_x_next = np.concatenate((opt_x_next, extension), axis=0)

                init_guess = np.concatenate((opt_x_next,opt_u_next), axis=0)

            opt_x, opt_u = self.ocp.solve(self.obj.x, init_guess,self.solver)
            u_k = opt_u[0:ocp_steps_per_mpc_iter]

            if i == 0:
                plot_trajectories(opt_x.T, opt_u, self.obj.h, title="OCP Iteration 1")

            for u_i in u_k:

                # Applying first part of optimal input trajectory
                self.obj.x = self.obj.step(u_i)

                # Uncomment the following lines to simulate arduino communication
                # self.ser.send(u_i)
                # self.obj.x = self.ser.receive().reshape(self.obj.nx, 1)

                # Appending to MPC result
                self.res_x_mpc.append(self.obj.x)
                self.res_u_mpc.append(u_i)

        self.res_x_mpc = np.concatenate(self.res_x_mpc, axis=1)
        self.res_u_mpc = np.array(self.res_u_mpc)

        return self.res_x_mpc,self.res_u_mpc


class SimulatedInvertedPendulum:
    def __init__(self, g, l, b, m, h, x0):
        self.nx = 4  # Number of states
        self.nu = 1  # Number of inputs

        self.x_min = [-np.inf, -np.inf, -0.3, -np.inf]  # State lower limits
        self.x_max = [np.inf, np.inf, 0.3, np.inf]  # State upper limits

        self.u_min = [-23]  # Input lower bound
        self.u_max = [23]  # Input upper bound

        self.x = x0.reshape(self.nx, 1) # Initial state as a numpy array
        
        self.h = h # Discretisation rate
        
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
                       (3*self.g/(2*self.l))*np.sin(x[0]) - (3*self.b/(self.m*self.l**2))*x[1] + 3/(2*self.l)*np.cos(x[0])*u,
                       x[3],
                       u
                       )
        return xdot
       
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
        k2 = self.system(self.x+self.h/2*k1, u)
        k3 = self.system(self.x+self.h/2*k2, u)
        k4 = self.system(self.x+self.h*k3,  u)
        self.x = self.x + self.h/6 * (k1 + 2*k2 + 2*k3 + k4)

        return self.x
    
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


def main():
    # ser = SerialCom('COM8', 256000)

    x_0 = np.array([pi, 0, 0, 0])
    simulation = True

    if simulation:
        pend = SimulatedInvertedPendulum(g=9.81,
                                         l=0.41,
                                         b=0.016,
                                         m=0.3,
                                         h=1/125,
                                         x0=x_0)
        ocp = OCP(1.5, pend.h, pend, lqr)
        dt = 1 / 12.5

        contr = Controller(ocp, pend, dt)

        x_mpc, u_mpc = contr.mpc_sim(T_sim=10)

        plot_trajectories(x_mpc,u_mpc,pend.h,'MPC')

    plt.show()


if __name__ == "__main__":
    main()