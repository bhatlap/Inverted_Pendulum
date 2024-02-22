%% Setup 
clear all
clc

syms x xd th thd xdd 
syms m l g b 

%%

xx = [th; thd; x; xd];
uu = [xdd];

thdd = (3*cos(th)/(2*l))*xdd + (3*sin(th)/(2*l))*g - (3*b/(m*l^2))*thd;

f = [thd; thdd; xd; xdd];

Jx = jacobian(f, xx);
Ju = jacobian(f, uu);

A = subs(Jx, [th, thd], [0,0]);
B = subs(Ju, [th, thd], [0,0]);

% (b / m) must be equal to 0.0533 according to the existing linear model

mm = 0.3;        % mass of the pendulum [kg]      
gg = 9.80665;    % gravitational acceleration [m / s^2]
bb = 0.016;      % friction coefficient
ll = 0.4;        % length of the pendulum [m]

A_assumed = subs(A, [m, g, b, l], [mm, gg, bb, ll]);
B_assumed = subs(B, [m, g, b, l], [mm, gg, bb, ll]);


