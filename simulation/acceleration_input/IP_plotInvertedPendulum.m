
function IP_plotInvertedPendulum(t, XX, UU)
    line_width = 1.5;
    
    figure;
    sgtitle('States') 
    
    subplot(2,2,1)
    plot(t, XX(1,:),'Linewidth', line_width)
    xlabel('time (s)')
    ylabel('th (rad)')
    
    subplot(2,2,2)
    plot(t, XX(2,:),'Linewidth', line_width)
    xlabel('time (s)')
    ylabel('thd (rad/s)')
    
    subplot(2,2,3)
    plot(t, XX(3,:),'Linewidth', line_width)
    xlabel('time (s)')
    ylabel('x (m)')
    
    subplot(2,2,4)
    plot(t, XX(4,:),'Linewidth', line_width)
    xlabel('time (s)')
    ylabel('xd (m/s)')

    figure;
    sgtitle('Inputs')
    stairs(t, UU, 'Linewidth', line_width)
    xlabel('time (s)')
    ylabel('xdd (m/s^2)')
    
end
    