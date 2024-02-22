%% NMPC -- Theory and Applications
% Course at Karlsruhe Institute of Technology
% Problem Set #3 Exercise 1
%
%
% Author(s):    https://github.com/casadi/optistack/blob/master/examples/rk4.m
% Date:         April 10, 2017
%% Runge-Kutte 4 Integrator
function xf = rk4(ode,h,t,x,u)
  k1 = ode(t,x,u);
  k2 = ode(t,x+h/2*k1,u);
  k3 = ode(t,x+h/2*k2,u);
  k4 = ode(t,x+h*k3,  u);
  xf = x + h/6 * (k1 + 2*k2 + 2*k3 + k4);
end