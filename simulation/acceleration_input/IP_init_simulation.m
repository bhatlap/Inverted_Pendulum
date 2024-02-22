%% Setup

clear all
clc

%% Simulation Parameters 

t0 = 0;                        % start time
t1 = 10;                       % stop time
    
T = 0.008;                     % sampling time

t = t0:T:t1;                   % time horizon

sim_len = size(t,2);           % simulation length                 

%% Model Parameters

params.m = 0.3;        % mass of the pendulum [kg]      
params.g = 9.80665;    % gravitational acceleration [m / s^2]
params.b = 0.016;      % friction coefficient
params.l = 0.4;        % length of the pendulum [m]

%% Initial Conditions

%th_init  = 0.1745;       % initial pendulum angle (5 degrees = 0.1745)
th_init  = pi;           % initial pendulum angle

thd_init = 0;            % initial pendulum velocity 
x_init   = 0;            % initial cart position
xd_init  = 0;            % initial cart velocity
             
init_x = [th_init; thd_init; x_init; xd_init]; % initial states
