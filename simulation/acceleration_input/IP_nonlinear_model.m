
function xxd = IP_nonlinear_model(xx, uu, params)

    m = params.m;
    g = params.g;
    l = params.l;
    b = params.b;
    
    th  = xx(1);         % pendulum's angle [rad]
    thd = xx(2);         % pendulum's angular velocity [rad/s]
    x   = xx(3);         % cart's position [m]
    xd  = xx(4);         % cart's velocity [m/s]
    
    xdd = uu;            % acceleration of the cart [m/s^2]

    thdd = (3*cos(th)/(2*l))*xdd + (3*sin(th)/(2*l))*g - (3*b/(m*l^2))*thd;

    xxd = [thd; thdd; xd; xdd];
end