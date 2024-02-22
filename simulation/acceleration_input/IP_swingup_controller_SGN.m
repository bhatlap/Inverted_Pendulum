
function uu = IP_swingup_controller_SGN(xx, x_star, params)

    m = params.m;
    g = params.g;
    l = params.l;
    
    th  = xx(1);
    thd = xx(2);
    x = xx(3);
    
    J = (m*l^2)/3;
    
    E_star = m*g*l;
    E_pendulum = 0.5*J*(thd^2) + m*g*(l/2)*(cos(th) + 1);
    
    kE = 3;
    kC = 20;
   
    % error calculations
    EE = E_star - E_pendulum;
    PE = x_star - x;
    
    if EE*thd*cos(th) == 0
        uu = kE + kC*PE;
    else
        uu = kE*sign(EE*thd*cos(th)) + kC*PE;
    end
    
    if uu > 25
        uu = 25;
    elseif uu < -25
        uu = -25;
    end
    
end

