%% Setup

IP_init_simulation;
 
%% Simulation

XX = zeros(size(init_x,1), size(t,2));
UU = zeros(1, size(t,2));

XX(:,1) = init_x;

for i = 1:size(t,2)-1
   
     xx = XX(:,i);
     
     uu = IP_stabilization_controller(xx, T, params);
     UU(:, i) = uu; 
     
     xf = rk4(@(x, u) IP_nonlinear_model(x, u, params), T, xx, uu);
     XX(:, i+1) = xf;
     
end

IP_plotInvertedPendulum(t, XX, UU)

