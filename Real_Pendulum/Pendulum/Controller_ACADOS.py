
from acados_template import AcadosOcp, AcadosOcpSolver, AcadosSimSolver, AcadosSim, AcadosModel
from Pendulum.Model import pend_model, RealInvertedPendulum, SimulatedInvertedPendulum
from Pendulum.LQR import StabilizationLQR
from casadi import *
import matplotlib.pyplot as plt
import pandas as pd
import scipy.linalg as lin


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


class OCP:
    """
    Class to interact with and define the Optimal Control Problem.
    """

    def __init__(self, xr=np.array([0, 0, 0, 0]), horizon=1, f=20,port='/dev/ttyUSB0', baud_rate=115200):

        self.ocp = AcadosOcp()  # OCP Object to formulate OCP
        self.ocp.model = pend_model()  # Pendulum model
        self.ocp.code_export_directory = "export_dir"

        self.nx = self.ocp.model.x.size()[0]
        self.nu = self.ocp. model.u.size()[0]
        self.ny = self.nx + self.nu

        self.xr = xr  # Set-point/ Goal
        self.Q = np.diag([1e1, 1e-4, 1e4, 1e-2])
        self.R = 1e1 * np.eye(self.nu)

        self.Tf = horizon  # Prediction horizon in seconds
        self.N = f  # Steps per prediction horizon (discretization frequency)
        self.counter = 0

        self.port = port
        self.baud_rate= baud_rate

    def build(self, simulation):
        """
        This function builds an OCP of the form:

            min(x,u) Σ(xT Q x + uT R u) + xT Qt x
                s.t
                    x_dot = model
                    xlb <= x <= xub
                    ulb <= u <= uub

                x : [theta dtheta/dt p dp/dt]

        :return: ocp_solver : The solver object to be used to solve the OCP
        """

        if not simulation:
            port = self.port
            baud_rate = self.baud_rate
            pend = RealInvertedPendulum(g=9.81, l=0.41, b=0.016, m=0.3, h=self.Tf/self.N, port=port, baud_rate=baud_rate)
            lqr = StabilizationLQR(pend)
            Q_ter, _ = lqr.lqr()

        else:
            x_init = np.array([np.pi, 0, 0, 0]).reshape(self.nx,1)
            pend = SimulatedInvertedPendulum(g=9.81, l=0.41, b=0.016, m=0.3, h=self.Tf/self.N, x0=x_init)
            lqr = StabilizationLQR(pend)
            Q_ter, _ = lqr.lqr()

        """ OCP MODELING """

        """ Cost Function """
        self.ocp.cost.cost_type = 'LINEAR_LS'
        self.ocp.cost.cost_type_e = 'LINEAR_LS'
        self.ocp.cost.W = lin.block_diag(self.Q, self.R)
        self.ocp.cost.W_e = Q_ter
        self.ocp.cost.Vx = np.zeros((self.ny, self.nx))
        self.ocp.cost.Vx[:self.nx, :self.nx] = np.eye(self.nx)
        Vu = np.zeros((self.ny, self.nu))
        Vu[4, 0] = 1.0
        self.ocp.cost.Vu = Vu
        self.ocp.cost.Vx_e = np.eye(self.nx)
        self.ocp.cost.yref = np.append(self.xr, 0.)
        self.ocp.cost.yref_e = self.xr

        """ State Constraints """
        self.ocp.constraints.constr_type = 'BGH'
        self.ocp.constraints.lbx = np.array(pend.x_min[2:])  # [-0.28,-1.8]
        self.ocp.constraints.ubx = np.array(pend.x_max[2:])  # [0.19,1.8]
        self.ocp.constraints.idxbx = np.array([2,3])  # Constraint applied to 3rd state (cart position)

        """Input Constraints"""
        self.ocp.constraints.lbu = np.array(pend.u_min)
        self.ocp.constraints.ubu = np.array(pend.u_max)
        self.ocp.constraints.idxbu = np.array([0])
        self.ocp.constraints.x0 = np.array([np.pi, 0.0, 0.0, 0.0])

        """Solver Options"""
        self.ocp.dims.N = self.N
        self.ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
        self.ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
        self.ocp.solver_options.integrator_type = 'IRK'  # System model integrator
        self.ocp.solver_options.print_level = 0
        self.ocp.solver_options.nlp_solver_type = 'SQP_RTI'
        self.ocp.solver_options.qp_solver_cond_N = self.N  # Stages per Prediction horizon
        self.ocp.solver_options.tol = 1e-2  # Tolerance
        self.ocp.solver_options.nlp_solver_max_iter = 100
        self.ocp.solver_options.tf = self.Tf  # Prediction horizon
        self.ocp.solver_options.levenberg_marquardt = 10.0

        """ Generate .json file to build solver """
        AcadosOcpSolver.generate(self.ocp, json_file='acados_MPC.json')
        AcadosOcpSolver.build(self.ocp.code_export_directory, with_cython=True)

    def solve(self, solver):
        """
        :param x: Current state of the system (i.e. next initial state)
        :param init_x: Initial guess to the solver
        :param init_u: Initial guess to the solver
        :param solver: Solver object
        :return: opt_x, opt_u : optimal x and u trajectories

        This function solves the OCP and returns an optimal state and input trajectory
        """
        """ Solve the OCP """
        result = solver.solve()

        """ Extract Solution """
        N = self.N
        nx = self.nx

        opt_x = np.zeros([nx, N + 1])
        for i in range(N + 1):
            opt_x[:, i] = solver.get(i, 'x')

        opt_u = np.zeros([1, N])
        for i in range(N):
            opt_u[0, i] = solver.get(i, 'u')

        return opt_x.T, opt_u.T, result


class MPC:
    """
    Class to interact with the Model Predictive Controller
    """
    def __init__(self, pendulum_object, xr=np.array([0, 0, 0, 0]), Tf=1, freq=20, port='/dev/ttyUSB0', baud_rate=115200, simulation=False):

        self.ocp = OCP(xr=xr, horizon=Tf, f=freq, port= port, baud_rate=baud_rate)  # Object of class OCP
        self.obj = pendulum_object  # Object of class defining the system

        #self.ocp.build(simulation=simulation)  # Builds .json files to create Solver objects
                                                 # (Uncomment this line if you make changes to the OCP)

        self.ocp_solver = AcadosOcpSolver.create_cython_solver(json_file='acados_MPC.json')
        self.dt = self.ocp.Tf/self.ocp.N  # Sampling time of controller
        self.opt_x = None
        self.opt_u = None
        self.first_time = True
        self.counter = 0
        self.res = []

    def set_init_guess(self, init_x, init_u):
        """
        This function sets the initial guess for the solver to start optimization
        :param init_x: Initial guess for the states (dimensions: ((N+1)*nx,1))
        :param init_u: Initial guess for the inputs (dimensions: (N*nu,1))
        :return: sets the initial guess for the solver
        """

        for i in range(self.ocp.N + 1):
            self.ocp_solver.set(i, 'x', init_x[i, :])

        for i in range(self.ocp.N):
            self.ocp_solver.set(i, 'u', init_u[i])

    def set_init_state(self,x):
        """
        This function sets the initial state for the start of the state trajectory
        :param x: State to be set as the initial state
        :return: sets the state 'x' as the initial state for the state trajectory
        """
        x = x.reshape(self.ocp.nx,)

        self.ocp_solver.set(0, "lbx", x)
        self.ocp_solver.set(0, "ubx", x)

    def step(self, x):
        """
        This function calculates the next part of the input trajectory by solving the optimization problem

        :param x: Current state
        :return u_k: First part of the input trajectory
        """

        """ Generating initial guess for solver """
        if self.first_time:
            # First initial guess
            # init_x, init_u = (np.repeat(x, self.ocp.N + 1, axis=1).T, [0] * self.ocp.N)
            # init_u = np.array(init_u).reshape(self.ocp.N*self.ocp.nu,1)
            init_x = pd.read_csv('Pendulum/init_x.csv').to_numpy()[:, 1:]
            init_u = pd.read_csv('Pendulum/init_u.csv').to_numpy()[:, 1:]

        else:
            # Warm-starting the solver
            init_u = self.opt_u.flatten()[self.obj.nu:self.obj.nu * self.ocp.N]
            init_u = np.concatenate((init_u, self.opt_u[self.ocp.N - 1]), axis=0).reshape(self.ocp.N * self.obj.nu, 1)

            init_x = self.opt_x[1:self.obj.nx * (self.ocp.N + 1), :]
            extension = self.obj.simulate(self.opt_x[self.ocp.N], self.opt_u[self.ocp.N - 1]).full().reshape(1,self.ocp.nx)
            init_x = np.concatenate((init_x, extension), axis=0)

        """ Set initial guess and state for solver """
        self.set_init_guess(init_x=init_x, init_u=init_u)
        self.set_init_state(x)

        """ Solving OCP """
        self.opt_x, self.opt_u, result = self.ocp.solve(self.ocp_solver)
        self.counter += 1
        self.res.append(result)

        if self.first_time:
            plot_trajectories(self.opt_x.T, self.opt_u, self.dt, 'Iteration 1')

        self.first_time = False

        u_k = self.opt_u[0]
        x_k = self.opt_x[1, :]  # Next predicted state

        return u_k, x_k, self.opt_u, self.opt_x






