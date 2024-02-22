from casadi import *
import numpy as np
import scipy.linalg as lin
from scipy import signal


class StabilizationLQR:
    def __init__(self, pend):
        self.pend = pend
        self.Q = np.diag([10, 0.1, 10, 0.1])
        self.R = 1 * np.eye(self.pend.nu)

        self.P = 0
        self.B = 0
        self.K = 0

    def lqr(self):

        x = SX.sym('x', self.pend.nx, 1)
        u = SX.sym('u', self.pend.nu, 1)
        xdot = self.pend.system(x, u)
        nx = self.pend.nx
        nu = self.pend.nu

        A_fun = Function('Linearized_A', [x, u], [jacobian(xdot, x)])
        B_fun = Function('Linearized_B', [x, u], [jacobian(xdot, u)])

        xs = np.zeros([nx, 1])
        us = np.zeros([nu, 1])

        A = A_fun(xs, us)
        B = B_fun(xs, us)
        C = np.eye(nx)
        D = np.zeros([nx, nu])

        sys = signal.StateSpace(A, B, C, D)
        sys_d = sys.to_discrete(dt=self.pend.h)
        A_d = sys_d.A
        B_d = sys_d.B

        P = lin.solve_discrete_are(A_d, B_d, self.Q, self.R)

        return P, B_d

    def step(self, x):
        self. P, self.B = self.lqr()
        self.K = 1/self.R * self.B.T @ self.P
        u = - np.dot(self.K, x)

        if u < -23:
            u = -23
        elif u > 23:
            u = 23

        return u

