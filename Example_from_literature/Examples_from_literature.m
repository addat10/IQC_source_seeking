close all
clear
clc
%% Vehicle Dynamics
% Example from Scherer, Weiland LMI notes (not a vehicle)
% a=1.2;
% G_veh=ss(tf([1,-a,0],[1,1,2,1]));

% Example motivated from Freeman(not a vehicle)
% G_veh=d2c(-tf([2 -1],[20 -10 10],1));

% Example from Freeman(not a vehicle)
% G_veh=d2c(-tf([2 -1],[20 -10 10],0.1));

% Practical example of a generic vehicle
% G_veh=tf([1],[1,1,0]);

% Example from Chen and Wen(not a vehicle)
G_veh=tf([3,3],[1,1,25,0,0]);

% Example modified from Chen and Wen(not a vehicle)
% G_veh=tf([3,3],[1,1,25,0]);

%% Verify exponential stability: Analysis
addpath('.\analysis_scripts')
m=1;
L=1.305; 

cvx_tol=1e-6;
bisect_tol=1e-2;
alpha_lims=[1e-6,10]; 
cond_tol=100000000;

% Circle criterion
alpha=1e-6;
%[status,P]=verify_exp_stab_CC(G_veh,alpha,m,L,cond_tol,cvx_tol);

% Zames Falb Multipliers
%alpha=0.2;
rho=-1;
psi_order=1;
[status,P]=verify_exp_stab_ZF_basis(G_veh,alpha,m,L,rho,psi_order,cond_tol,cvx_tol);
