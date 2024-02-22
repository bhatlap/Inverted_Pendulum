function IP_plotEnergy(t, XX, params)

    m = params.m;          % mass of the pendulum [kg]
    g = params.g;          % gravitational acceleration [m/s^2]
    l = params.l/2;        % half-length of the pendulum [m]
    
    E_star = (2*m*g*l)*ones(size(t));
    E_pendulum = zeros(size(t));
    
    for i = 1:size(t,2)
       th  = XX(1, i);
       thd = XX(2, i);
       
       E_pendulum(i) = 0.5*m*(l^2)*(thd^2) + m*g*l*(cos(th) + 1);
    end
    
    line_width = 1.5;
    
    figure;
    title('Pendulum''s Energy vs Reference Energy') 
    
    hold on
    plot(t, E_star, 'k', 'Linewidth', line_width)
    plot(t, E_pendulum, 'r', 'Linewidth', line_width)
    hold off
    xlabel('time (s)')
    ylabel('Energy')
    legend('Reference Energy', 'Pendulum Energy')
    
end