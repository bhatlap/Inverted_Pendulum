function [x,u] = MPCloop(U,X,X0,ocp,myfun,dt,tf,xinit,XX,UU)
    t           =   0:dt:tf;
    N           =   length(t)-1;
    nx          =   length(xinit);
    [nu,Nmpc]   =   size(U);
    x           =   zeros(nx,N+1);
    u           =   zeros(nu,N);
    x(:,1)      =   xinit;
%     Xsol        =   repmat(xfin,1,Nmpc+1);
%     Usol        =   repmat(uinit,1,Nmpc);
    Xsol        =   [XX(:,1) XX(:,1:Nmpc)];
    Usol        =   [UU(:,1) UU(:,1:Nmpc-1)];
    
    for k = 1:N
        % set current state
        ocp.set_value(X0, x(:,k)); 
        % initialize OCP
        ocp.set_initial(X,[Xsol(:,2:end) Xsol(:,end)]);
        ocp.set_initial(U,[Usol(:,2:end) Usol(:,end)]);
        % resolve OCP
        sol = ocp.solve();
        Xsol    =   sol.value(X);
        Usol    =   sol.value(U);
        % get first element
        u(:,k)   =   Usol(:,1);
        % apply to system using rk4
       x(:,k+1) =   rk4(@(t,x,u)myfun(t,x,u),dt,t(k),x(:,k),u(:,k));
        % apply to system using Matlab's ODE solvers
 %       [~,xt]     =   ode45(@(t,x,u)myfun(t,x,u),[t(k), t(k)+dt],x(:,k),[],u(:,k));
%         [~,xt]     =   ode15s(@(t,x,u)CSTR_ode(t,x,u),[tsim(k), tsim(k)+dt],xsim(:,k),[],usim(:,k));
  %      x(:,k+1) =   xt(end,:)';
    end
end