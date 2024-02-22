clear all;
close all;
clc;

t0   =  0;  % initial time
tf   =  10;  % end time
dt    =  0.008;
N    =  (tf-t0)/dt; % horizon length
Nmpc    =   (2-t0)/dt;
Tmpc    =   Nmpc*dt;
tmpc    =   0:dt:Tmpc;

t    =  t0:dt:tf;
x0   =  [pi; 0; 0; 0];
 m = 0.3;
 l = 0.4;
 g = 9.81;
nx = 4; 
nu = 1; 

load initialguess.mat;

Q   =   diag([1 1 1 1]);
R   =   eye(nu);
Qt  = diag([100 100 100 100]);
% constraints
umin    =   -25;
umax    =   25;

import casadi.*
ocp = casadi.Opti();
X       =   ocp.variable(nx,Nmpc+1);
U       =   ocp.variable(nu,Nmpc);
X0      =   ocp.parameter(nx);
J       =   0;
Xbar    =   [0; 0; 0; 0];%[1e-6; 1e-6; 1e-6; 1e-6];
for i=1:Nmpc
    X_next = rk4(@(t,x,u)ODE(t,x,u),dt,tmpc(i),X(:,i),U(:,i));
    %X_next = X(:,i) + dt*ODE(t(i),X(:,i),U(:,i));
    
    % gap closing constraint
    ocp.subject_to(X_next==X(:,i+1));

    ocp.subject_to(-0.42<=X(3,i)<=0.42);
    ocp.subject_to(umin<=U(:,i)<=umax);
    % cost function
    J=J+(X(:,i)-Xbar)'*Q*(X(:,i)-Xbar)+ U(:,i)'*R*U(:,i);%X(:,i)'*Q*X(:,i)

end
J=J+ (X(:,end)-Xbar)'*Qt*(X(:,end)-Xbar);
%ocp.subject_to(umin<=U<=umax);
ocp.set_value(X0, x0);
% state constraint for the last time step
ocp.subject_to(X(:,end)==rk4(@(t,x,u)ODE(t,x,u),dt,tmpc(end),X(:,end),U(:,end)));
%ocp.subject_to(X(:,end)==X(:,N) + dt*ODE(t(N),X(:,N),U(:,N)));

% initial condition
ocp.subject_to(X(:,1)==X0);
% initial guess
ocp.set_initial(X,repmat(Xbar,1,Nmpc+1));
%ocp.set_initial(U,repmat(0,1,N));
% set solver
ocp.solver('ipopt');

% set objective
ocp.minimize(J);

% solve ocp
%sol = ocp.solve();
[Xsol,Usol] = MPCloop(U,X,X0,ocp,@(t,x,u)ODE(t,x,u),dt,tf,x0,XX,UU);

% get solution
 %Xsol    =   sol.value(X);
 %Usol    =   sol.value(U);

figure
subplot(3,1,1)
plot(t, Xsol(1,:),'Linewidth', 2)
xlabel('t')
ylabel('x_1')
grid on

subplot(3,1,2)
stairs(t(1:N), Usol(:),'Linewidth', 2)
xlabel('t')
ylabel('u')

subplot(3,1,3)
stairs(t, Xsol(3,:),'Linewidth', 2)
xlabel('t')
ylabel('x3')

% figure
% plot(t, Xsol(1,:),'Linewidth', 2)
% xlabel('t')
% ylabel('x_1')
% grid on

%axis([0, 1.4, -20, 20]);
% grid on
% xlim([t_plot(1), t_plot(end)]);

% figure
% for i =1:N
%     x2=Xsol(3,i)+(l*cos(Xsol(1,i)));
%     y2=Xsol(3,i)+(l*sin(Xsol(1,i)));
%     plot([0 y2],[Xsol(3,i) x2])
%  
%     pause(0.01)
% end
function xf = eulerf(f,h,t,x,u)
    xf = x+h*f(t,x,u);
end
function f=ODE(t,x,u)
     m = 0.3;
    l = 0.4;
    g = 9.81;
    b = 0.016;
    f1 = x(2);
    f2 = (3*g)/(2*l)*sin(x(1))-(3*b)/(m*l^2)*x(2)+(3)/(2*l)*cos(x(1))*u;
    f3 = x(4);
    f4 = u;
    f=[f1;f2;f3;f4];
end