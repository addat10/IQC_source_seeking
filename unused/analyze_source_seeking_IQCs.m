clear
clc
%% System Dynamics
G_veh=define_vehicle_model();
%% Verify exponential stability: Analysis
m=1;
L=2;
[Psi_GI,M]=define_ZF_multiplier(m,L,G_veh);
%alpha=0.3;
cvx_tol=1e-3;
bisect_tol=1e-2;
%alpha_lims=[0,10];
alpha_lims=[0.1,10]; % alpha_best=0.2744
[alpha_best,~]=bisection_exponent(Psi_GI,M,alpha_lims,cvx_tol,bisect_tol);
[status,P]=verify_exp_stab(Psi_GI,M,alpha_best,cvx_tol*10);
%% Numerically simulate the dynamics
y_min=1*(-1+2*rand());
k=2;
grad_field=@(y) k*(y-y_min);
sim_time=100;
dt=0.01;
time_steps=sim_time/dt;
pos_ic=10*(-1+2*rand());
x_ic=[pos_ic;0];
[trajs]= simulate_source_seek(G_veh,x_ic,grad_field,time_steps,dt);
%% Theoretical upper bound based on LMIs
time=dt*(1:time_steps);
e_ub=exp(-alpha_best*time)*cond(P)*norm(x_ic-[y_min;0])^2;
%% Asymptotic Lyapunov Exponent for known quadratic fields
A_cl=G_veh.A-G_veh.B*k*G_veh.C;
% Compute the known theoretical lower and upper bounds

% Lower bound: a factor 2 appears as we are estimating norm squares.
alpha_lb_eig_lr= -2*(real(eigs(A_cl, 1, 'lr'))); 
e_lb_eig_lr=exp(-alpha_lb_eig_lr*time)*norm(x_ic-[y_min;0])^2;

% Upper bound: a factor 2 appears as we are estimating norm squares.
sym_A_cl=0.5*(A_cl'+A_cl);
alpha_ub_eig_sym_A=-2*(eigs(sym_A_cl, 1, 'lr'));
e_ub_eig_sym_A=exp(-alpha_ub_eig_sym_A*time)*norm(x_ic-[y_min;0])^2;
%% Plot
figure()
plot(trajs.x(1,:))
%% Compare the obtained numerical decay with the theoretical decay
e=trajs.x(:,:)-[y_min;0];
e_norms=(ones(1,2)*(e.*e));

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

