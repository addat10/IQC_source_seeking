function []=animate_formation_source_seeking(traj,time_steps,L,X,Y,Z)
    % This function animates the trajectories of the agents    
    n=size(L,1);
    dim=size(traj,1)/n;
    figure()
    switch dim        
        case 2
            for time=1:time_steps
                if rem(time,100)==0
                    contour(X,Y,Z,20)
                    hold on
                    for i=1:n
                        xy=[(i-1)*2+1;(i-1)*2+2];
                        plot(traj(xy(1),time),traj(xy(2),time),'x')
                    end
                    hold off
                    pause(0.1)
                end
            end
                  
    end
end