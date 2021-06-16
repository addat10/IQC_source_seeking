clear
clc
%% System Dynamics
G=define_vehicle_model();
m=1;
L=2;
[Psi_GI,M]=define_ZF_multiplier(m,L,G)
%% Verify exponential stability
alpha=0.3;
cvx_tol=1e-3;
bisect_tol=1e-2;
alpha_lims=[0,10];
alpha_best=bisection_exponent(Psi_GI,M,alpha_lims,cvx_tol,bisect_tol)
status=verify_exp_stab(Psi_GI,M,alpha_best,cvx_tol*10)