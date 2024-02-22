%% Setup

IP_init_simulation;

%% Trajectories

XX = zeros(size(init_x,1), sim_len);
XX(:,1) = init_x;

UU = zeros(1, size(t,2));

x_star = 0;
iEE = 0;

%% Simulation

for i = 1:size(t,2)-1
    
    xx = XX(:,i);
    
    uu = IP_swingup_controller_SGN(xx, x_star, params);
    %uu = IP_swingup_controller_new(xx, x_star, params);
    %[uu, iEE] = IP_swingup_controller(xx, iEE, T, x_star, params);
    UU(i) = uu;
    
    xf = rk4(@(x, u) IP_nonlinear_model(x, u, params), T, xx, uu);
    XX(:, i+1) = xf;
    
end

%% Plot

IP_plotInvertedPendulum(t, XX, UU)
IP_plotEnergy(t, XX, params)

