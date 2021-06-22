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
alpha_lims=[0,10];
[alpha_best,~]=bisection_exponent(Psi_GI,M,alpha_lims,cvx_tol,bisect_tol);
[status,P]=verify_exp_stab(Psi_GI,M,alpha_best,cvx_tol*10);
%% Numerically simulate the dynamics
y_min=10*(-1+2*rand());
grad_field=@(y) 1*(y-y_min);
sim_time=100;
dt=0.01;
time_steps=sim_time/dt;
pos_ic=10*(-1+2*rand());
x_ic=[pos_ic;0];
[trajs]= simulate_source_seek(G_veh,x_ic,grad_field,time_steps,dt);
%% Theoretical upper bound
time=dt*(1:time_steps);
e_ub=cond(P)*norm(x_ic-[y_min;0])*exp(-alpha_best*time);
%% Plot
figure()
plot(trajs.x(1,:))
%% Compare the obtained numerical decay with the theoretical decay
e=trajs.x(:,:)-[y_min;0];
e_norms=sqrt(ones(1,2)*(e.*e));
figure()
plot(time,log(e_norms))
hold on
plot(time,log(e_ub))
ylim([-10,2])
legend('e norm','theoretical bound')
xlabel('time')
ylabel('ln(pos error)')

