function [trajs]= simulate_source_seek(G_veh,x_ic,grad_field,time_steps,dt)
% This function simulates the feedback loop (G,grad)
    [A,B,C,D]=ssdata(G_veh);
    nx=size(A,1);
    [ny,nu]=size(D);
    x=zeros(nx,time_steps);
    y=zeros(ny,time_steps);
    u=zeros(ny,time_steps);    
    if D~=zeros(ny,nu)
        error('Vehicle model not strictly proper. (Necessary to avoid algebraic loop)')
    end
    % Initial condition
    x(:,1)=x_ic;    
    y(:,1)=C*x_ic;
    for i=1:(time_steps-1)
        u(:,i)=-grad_field(y(:,i));
        x(:,i+1)=(eye(nx)+dt*A)*x(:,i)+dt*B*u(:,i);
        y(:,i+1)=C*x(:,i);
    end
    trajs=struct;
    trajs.x=x;
    trajs.y=y;
    trajs.u=u;
end