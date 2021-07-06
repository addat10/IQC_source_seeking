function [trajs]= simulate_formation_source_seek(G_veh,x_ic,form_ref,grad_field,Lap,time_steps,dt)
% This function simulates the feedback loop (G,grad)
    [A_veh,B_veh,C_veh,D_veh]=ssdata(G_veh);
    
    % Build the full system
    n=size(Lap,1);
    A=kron(eye(n),A_veh);
    B=kron(eye(n),B_veh);
    C=kron(eye(n),C_veh);
    D=kron(eye(n),D_veh);
    
    % Get sizes of signals
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
    dim=ny/n;
    for i=1:(time_steps-1)
        u(:,i)=kron(Lap,eye(dim))*(form_ref-y(:,i));
        u(1:dim,i)=u(1:dim,i)-grad_field(y(1:dim,i));
        x(:,i+1)=(eye(nx)+dt*A)*x(:,i)+dt*B*u(:,i);
        y(:,i+1)=C*x(:,i);
    end
    trajs=struct;
    trajs.x=x;
    trajs.y=y;
    trajs.u=u;
end