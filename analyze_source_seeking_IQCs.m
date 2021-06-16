clear
clc
%% System Dynamics
G=define_vehicle_model();
m=1;
L=2;
[Psi_GI,M]=define_ZF_multiplier(m,L,G);
%% Verify exponential stability
%alpha=0.3;
cvx_tol=1e-3;
bisect_tol=1e-2;
alpha_lims=[0,10];
[alpha_best,P]=bisection_exponent(Psi_GI,M,alpha_lims,cvx_tol,bisect_tol);
[status,P]=verify_exp_stab(Psi_GI,M,alpha_best,cvx_tol*10);
%% Simulate the dynamics
y_min=2;
grad_field=@(y) 2*(y-y_min);
time_steps=10000;
dt=0.01;
pos_ic=10*(-1+2*rand());
x_ic=[pos_ic;0];
[trajs]= simulate_source_seek(G,x_ic,grad_field,time_steps,dt);
%% Theoretical upper bound
time=dt*(1:time_steps);
e_ub=cond(P)*norm(x_ic-[y_min;0])*exp(-alpha_best*time);
%% Plot
figure()
plot(trajs.x(1,:))
%%
e=trajs.x(:,:)-[y_min;0];
e_norms=sqrt(ones(1,2)*(e.*e));
figure()
plot(e_norms)
hold on
plot(e_ub)
legend('e norm','theoretical bound')

