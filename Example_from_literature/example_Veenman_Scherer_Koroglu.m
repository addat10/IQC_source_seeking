close all
clear
clc
%% Plant Dynamics
% Example from Veenman, Scherer and Koroglu (not a vehicle)
% A=[-0.4,    -1;...
%     1,      0];
% B=[-0.2,    -1,     -0.25;...
%     0,      0,      0];
% C=[1,   0;...
%    0,   1;...
%    0,   0];
% D=zeros(3);D(3,2)=1;
% Example from Fetzer and Scherer page 3389 Example 4.9 (not a vehicle)
A=[-4,      -3,     0;...
    2,      0,      0;...
    -1,     -1,     -2];
B=-[0,   4,      1,      3;...
   2,   0,      3,      1;...
   1,   0,      3,      1];
C=[-0.1,    -0.2,      1;...
   -1,      -0.3,      0.1;...
   -0.2,    0.1,       1;...
   0.1,    -0.2,       0.2;];
D=zeros(4);
% 
G_veh=ss(A,B,C,D);

%% Verify exponential stability: Analysis
addpath('.\analysis_scripts')
m=0;
L=2.1/4; 

cvx_tol=1e-6;
bisect_tol=1e-2;
cond_tol=1000000000;
alpha=0;

% Circle criterion
%[status,P]=verify_exp_stab_CC_MIMO(G_veh,alpha,m,L,cond_tol,cvx_tol);

%Full Block Circle criterion
%[status,P]=verify_exp_stab_FBCC(G_veh,alpha,m,L,cond_tol,cvx_tol);


% Zames Falb Multipliers
rho=-1;
psi_order=1;
[status,P]=verify_exp_stab_ZF_basis_3(G_veh,alpha,m,L,rho,psi_order,cond_tol,cvx_tol);

% rho=-100;
% psi_order=1;
% [status,P]=verify_exp_stab_ZF_basis(G_veh,alpha,m,L,rho,psi_order,cond_tol,cvx_tol);