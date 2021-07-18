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
        c_damp=2;mass=1;step_size=1;
        dim=2;% spatial dimension (of positions and velocities)
        G_veh=define_G_mass_with_friction_wrapped(dim,c_damp,mass,step_size);
    case 2
        % Quadrotor dynamics
        addpath(genpath('.\vehicles\quadrotor'))
        dim=2;% spatial dimension (of positions and velocities)
        % Current implementation only supports dim=2 for quadrotors
        G_veh=define_G_quad_wrapped(dim);        
end
%% Network Properties
n=5; % Number of agents
link_prob=0.5; topo='rand';
A=gen_topology(n,link_prob,topo);
D=A*ones(n,1);
Lap=diag(D)-A;
eig_Ls=eig(Lap);
%% Field Properties
range=10;
y_min=range*(-1+2*rand(dim,1));
alpha_bar=1; % eig_Ls(2)/4 maximizes the lower bound on the grounded Lap
beta_bar=1.5;
switch(dim)
    case 1
        k=min(alpha_bar,beta_bar);
        grad_field=@(y) k*(y-y_min);
    case 2
        Hess=[alpha_bar,0;0,beta_bar];
        x = linspace(-2*range,2*range);
        y = linspace(-2*range,2*range);
        [X,Y] = meshgrid(x,y);
        Z = 1*(X-y_min(1)).^2+2*(Y-y_min(2)).^2;
        grad_field=@(y) Hess*(y-y_min);
end
%% Verify exponential stability: Analysis
addpath('.\analysis_scripts')
m=(alpha_bar/n)*(1-(2*alpha_bar/eig_Ls(2)));
L_per=Lap;
L_per(1,1)=L_per(1,1)+alpha_bar;
m_true=min(eig(L_per));
L=max(eig_Ls)+beta_bar;
cvx_tol=1e-6;
bisect_tol=1e-2;
alpha_lims=[1e-9,10]; 
cond_tol=100;

% Multiplier class
% Select a multiplier class from the following choices
% 1. Circle criterion
% 2. Full block circle criterion
% 3. Zames Falb multipliers
multiplier_flag=3;
[alpha_best,~]=bisection_exponent(G_veh,m,L,alpha_lims,cond_tol,cvx_tol,bisect_tol,multiplier_flag);
[status,P]=verify_exp_stab_ZF(G_veh,alpha_lims(1),m,L,cond_tol,cvx_tol);
%% Numerically simulate the dynamics
% Define the underlying field for dim=2
sim_time=100;
dt=0.001;
time_steps=sim_time/dt;
switch(Veh_mod)
    case 1
        % Mass with friction dynamics
        x_ic=5*(-1+2*rand(n*2*dim,1));
        thetas=linspace(0,2*pi,n);
        offsets=[cos(thetas(1:end-1));sin(thetas(1:end-1))];        
        form_ref=[zeros(dim,1);reshape(offsets,[(dim*(n-1)),1])];
        [trajs]= simulate_formation_source_seek(G_veh,x_ic,form_ref,grad_field,Lap,time_steps,dt); 
    case 2
        % Quadrotor dynamics
        x_ic=0;   
        % Need to implement
end
%% Theoretical upper bound based on LMIs
x_eqm=trajs.x(:,end);
time=dt*(1:time_steps);
e_ub_LMI=exp(-alpha_best*time)*cond(P)*norm(x_ic-x_eqm)^2;
%% Asymptotic Lyapunov Exponent for known quadratic fields
% Compute the known theoretical lower and upper bounds
DELTA=kron(L,eye(dim));
DELTA(1:dim,1:dim)=DELTA(1:dim,1:dim)+Hess;
switch(dim)
    case 1
        A_cl=G_veh.A-G_veh.B*k*G_veh.C;
    case 2
        A_cl=G_veh.A-G_veh.B*DELTA*G_veh.C;
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
        for i=1:n
            xy=[(i-1)*2+1;(i-1)*2+2];
            plot(trajs.y(xy(1),:),trajs.y(xy(2),:))
            hold on
            plot(trajs.y(xy(1),1),trajs.y(xy(2),1),'o')
            plot(trajs.y(xy(1),end),trajs.y(xy(2),end),'X')
        end        
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
plot(time,e_ub_LMI)
legend('e norm','LMI_bound')
xlabel('time')
ylabel('pos error')

% Plot trajectories on a logarithmic scale
figure()
plot(time,log(e_norms))
hold on
plot(time,log(e_ub_LMI))
plot(time,log(e_lb_eig_lr))
plot(time,log(e_ub_eig_sym_A))
ylim([-50,50])
legend('e norm','LMI','lb:2*Re(lambda max (A))','ub:eig(sym(A)')
xlabel('time')
ylabel('ln(pos error)')

%% Animate the trajectories
plot_steps=10000;
animate_formation_source_seeking(trajs.y,plot_steps,Lap,X,Y,Z)
