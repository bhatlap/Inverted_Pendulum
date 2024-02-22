
function u = IP_stabilization_controller(x, T, params)

    g = params.g;         % gravitational acceleration [m/s^2]
    l = params.l/2;       % half-length of the pendulum [m]

    % linear model of the inverted pendulum
    A = [0, 1, 0, 0;
         3*g/4/l, -1, 0, 0;
         0, 0, 0, 1;
         0, 0, 0, 0];
    B = [0; 3/4/l; 0; 1];
    C = [1 0 0 0; 0 0 1 0];
    D = [0; 0];

    % conversion to state-space representation
    sys_cont = ss(A,B,C,D);

    % conversion to discrete-time domain
    sys_disc = c2d(sys_cont, T);

    % declaration of weighting matrices
    Q = diag([1e3 1 100 1]);
    R = 10;

    % discrete LQR controller gain
    K = dlqr(sys_disc.A, sys_disc.B, Q, R);
    
    % applied input
    u = -K*x;
    
    % saturation of the input
    if u > 25
        u = 25;
    elseif u < -25
        u = -25; 
    end

end