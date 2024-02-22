from casadi import *
from acados_template import AcadosModel
import time
from Pendulum.SerialComm import SerialCom


def pend_model(l=0.41, m=0.3, b=0.016, g=9.81, model_name="OCL_Pendulum") -> AcadosModel:
    """

        Function referred from acados/examples/getting_started/pendulum_model.py
        :argument:
            m : Mass of Pendulum
            l : Length of pendulum
            b : Coefficient of angular friction
            g : Acceleration due to gravity

        :return:
            System Model: AcadosModel() object
    """



    theta = SX.sym('Theta')  # Angular displacement
    dtheta = SX.sym('dtheta')  # Angular velocity
    p = SX.sym('x')  # Linear displacement of cart
    dp = SX.sym('dp')  # Linear velocity of cart
    u = SX.sym('u')  # Input (Acceleration of cart)

    x = vertcat(theta, dtheta, p, dp)  # Vector of states
    u = vertcat(u)

    theta_dot = SX.sym('Theta_dot')
    dtheta_dot = SX.sym('dtheta_dot')
    p_dot = SX.sym('p_dot')
    dp_dot = SX.sym('dp_dot')

    xdot = vertcat(theta_dot, dtheta_dot, p_dot, dp_dot)
    f_expl = vertcat(x[1],
                     (3 * g / (2 * l)) * sin(x[0]) - (3 * b / (m * l ** 2)) * x[1] + 3 / (2 * l) * cos(x[0]) * u,
                     x[3],
                     u
                     )

    f_impl = xdot - f_expl

    model = AcadosModel()
    model.f_impl_expr = f_impl
    model.f_expl_expr = f_expl
    model.x = x
    model.xdot = xdot
    model.u = u
    model.name = model_name

    return model


class AbstractInvertedPendulum:
    def __init__(self, g, l, b, m, h):
        self.nx = 3  # Number of states
        self.nu = 2  # Number of inputs

        self.x_min = []  # State lower limits
        self.x_max = []  # State upper limits

        self.u_min = []  # Input lower bound
        self.u_max = []  # Input upper bound

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