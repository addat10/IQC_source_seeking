close all
clear
clc
%% System Dynamics
dim=2; % spatial dimension (of positions and velocities)
%G_veh=define_vehicle_model(dim);
G_veh=define_G_quad_wrapped(dim);
%% Verify exponential stability: Analysis
m=1;
L=2;
[Psi_GI,M]=define_ZF_multiplier(m,L,G_veh,dim);
%alpha=0.3;
cvx_tol=1e-6;
bisect_tol=1e-2;
%alpha_lims=[0,10];
alpha_lims=[0.0001,10]; % alpha_best=0.2744
[alpha_best,~]=bisection_exponent(Psi_GI,M,alpha_lims,cvx_tol,bisect_tol);
[status,P]=verify_exp_stab(Psi_GI,M,alpha_best,cvx_tol*10);
%% Numerically simulate the dynamics
% Define the underlying field
range=10;
y_min=range*(-1+2*rand(dim,1));
k=1;
x = linspace(-2*range,2*range);
y = linspace(-2*range,2*range);
[X,Y] = meshgrid(x,y);
Z = (X-y_min(1)).^2+(Y-y_min(2)).^2;
grad_field=@(y) k*(y-y_min);


sim_time=100;
dt=0.001;
time_steps=sim_time/dt;
%pos_ic=10*(-1+2*rand(dim,1));
%vel_ic=2*(-1+2*rand(dim,1));
%x_ic=[pos_ic;vel_ic];
%[trajs]= simulate_source_seek(G_veh,x_ic,grad_field,time_steps,dt);
[trajs]= simulate_source_seek_quad(G_veh,grad_field,time_steps,dt);

%% Theoretical upper bound based on LMIs
x_eqm=trajs.x(:,end);
time=dt*(1:time_steps);
x_ic=0;
e_ub=exp(-alpha_best*time)*cond(P)*norm(x_ic-x_eqm)^2;
%% Asymptotic Lyapunov Exponent for known quadratic fields
A_cl=G_veh.A-G_veh.B*k*G_veh.C;
% Compute the known theoretical lower and upper bounds

% Lower bound: a factor 2 appears as we are estimating norm squares.
alpha_lb_eig_lr= -2*max(real(eig(A_cl))); 
e_lb_eig_lr=exp(-alpha_lb_eig_lr*time)*norm(x_ic-x_eqm)^2;

% Upper bound: a factor 2 appears as we are estimating norm squares.
sym_A_cl=0.5*(A_cl'+A_cl);
alpha_ub_eig_sym_A=-2*max(eig(sym_A_cl));
e_ub_eig_sym_A=exp(-alpha_ub_eig_sym_A*time)*norm(x_ic-x_eqm)^2;
%% Plot
figure()
switch dim
    case 1
        plot(time,trajs.x(1,:)) 
        legend('Trajectory')
    case 2
        plot(trajs.y(1,:),trajs.y(2,:))
        hold on
        plot(trajs.y(1,1),trajs.y(2,1),'o')
        plot(trajs.y(1,end),trajs.y(2,end),'X')
        contour(X,Y,Z,20)
        legend('Trajectory','y(0)','y(end)')        
end



%% Compare the obtained numerical decay with the theoretical decay
e=trajs.x(:,:)-x_eqm;
e_norms=sum(e.^2);

% Plot trajectories on a linear scale
figure()
plot(time,e_norms)
hold on
plot(time,e_ub)

legend('e norm','theoretical bound with LMIs')
xlabel('time')
ylabel('pos error')

% Plot trajectories on a logarithmic scale
figure()
plot(time,log(e_norms))
hold on
plot(time,log(e_ub))
plot(time,log(e_lb_eig_lr))
plot(time,log(e_ub_eig_sym_A))
ylim([-10,10])
legend('e norm','theoretical bound with LMIs','lb:2*Re(lambda max (A))','ub:eig(sym(A)')
xlabel('time')
ylabel('ln(pos error)')

