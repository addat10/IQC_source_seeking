close all
clear
clc
%% Vehicle Dynamics
% Select a Vehicle model from the following choices
% 1. Mass with friction
% 2. Linearized Quadrotor
Veh_mod=1;
switch(Veh_mod)
    case 1
        % Mass with friction dynamics
        addpath(genpath('.\vehicles\mass_with_friction'))
        dim=1;% spatial dimension (of positions and velocities)
        G_veh=define_G_mass_with_friction_wrapped(dim);
    case 2
        % Quadrotor dynamics
        addpath(genpath('.\vehicles\quadrotor'))
        dim=2;% spatial dimension (of positions and velocities)
        % Current implementation only supports dim=2 for quadrotors
        G_veh=define_G_quad_wrapped(dim);        
end
%% Verify exponential stability: Analysis
addpath('.\analysis_scripts')
m=1;
L=2;
[Psi_GI,M]=define_ZF_multiplier(m,L,G_veh,dim);
%alpha=0.3;
cvx_tol=1e-6;
bisect_tol=1e-2;
%alpha_lims=[0,10];
alpha_lims=[0.0001,10]; % alpha_best=0.2744

% with usual Circle Criterion
[alpha_best_CC,~]=bisection_exponent(Psi_GI,M,alpha_lims,cvx_tol,bisect_tol);
alpha_best=alpha_best_CC;
%when using hinf desing, alpha_best=1e-4;
[status,P]=verify_exp_stab(Psi_GI,M,alpha_best,cvx_tol*10);

if dim==2
    % With full block circle criterion
    [alpha_best_FBCC,~]=bisection_exponent_FBCC(Psi_GI,m,L,dim,alpha_lims,cvx_tol,bisect_tol);
    [status,P]=verify_exp_stab_FBCC(Psi_GI,alpha_best,m,L,dim,cvx_tol);
    alpha_best=max(alpha_best,alpha_best_FBCC);
end
%%
% With Zames Falb Multiplier
[alpha_best_ZF,~]=bisection_exponent_ZF(G_veh,m,L,dim,alpha_lims,cvx_tol,bisect_tol);
alpha_best=max(alpha_best,alpha_best_ZF);
[status,P]=verify_exp_stab_ZF(G_veh,alpha_best,m,L,dim,cvx_tol);

%% Numerically simulate the dynamics
% Define the underlying field for dim=2
range=10;
y_min=range*(-1+2*rand(dim,1));
switch(dim)
    case 1
        k=min(m,L);
        grad_field=@(y) k*(y-y_min);
    case 2
        x = linspace(-2*range,2*range);
        y = linspace(-2*range,2*range);
        [X,Y] = meshgrid(x,y);
        Z = 1*(X-y_min(1)).^2+2*(Y-y_min(2)).^2;
        grad_field=@(y) [m,0;0,L]*(y-y_min);
end
sim_time=100;
dt=0.001;
time_steps=sim_time/dt;
switch(Veh_mod)
    case 1
        % Mass with friction dynamics
        pos_ic=10*(-1+2*rand(dim,1));
        vel_ic=2*(-1+2*rand(dim,1));
        x_ic=[pos_ic;vel_ic];
        [trajs]= simulate_source_seek(G_veh,x_ic,grad_field,time_steps,dt);        
    case 2
        % Quadrotor dynamics
        x_ic=0;
        [trajs]= simulate_source_seek_quad(G_veh,grad_field,time_steps,dt);
end
%% Theoretical upper bound based on LMIs
x_eqm=trajs.x(:,end);
time=dt*(1:time_steps);
e_ub=exp(-alpha_best*time)*cond(P)*norm(x_ic-x_eqm)^2;
%% Asymptotic Lyapunov Exponent for known quadratic fields
% Compute the known theoretical lower and upper bounds
switch(dim)
    case 1
        A_cl=G_veh.A-G_veh.B*k*G_veh.C;
    case 2
        A_cl=G_veh.A-G_veh.B*[m,0;0,L]*G_veh.C;
end
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

