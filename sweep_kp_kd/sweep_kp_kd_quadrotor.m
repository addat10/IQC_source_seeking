% This script is used to reproduce data for Example 7_20 from Scherer and
% Weiland's LMI notes. It is an exercise with an odd sector bounded
% non-linearity
close all
clear
clc
addpath('..\analysis_scripts')
addpath(genpath('..\vehicles\quadrotor'))
save_data=1;
save_path='.\data_quadrotor\';
%% Setup optimization
kp=[10,20,50,70,100];
kd=[10:20:1000];
n_p=length(kp);
n_d=length(kd);

% Sector bounds
m=1; % lower bound on the sector
L=10;  % Upper bound on the sector

% Optimization properties
cvx_tol=1e-3;
bisect_tol=1e-3;
cond_tol=100000000;

% Multiplier class
% Select a multiplier class from the following choices
% 1. Circle criterion
% 2. Full block circle criterion
% 3. Zames Falb multipliers
% 6. Zames Falb multipliers with basis
multiplier_flag=1;
save_path=[save_path,'mult_flag_',num2str(multiplier_flag)];
alpha_lims=[0,10]; % Initial range for the bisection algorithm
alpha_best=zeros(n_p,n_d);
%% Run the bisection algorithm over a 2D grid (a X L)
for i=1:n_p
    for j=1:n_d        
        kp_curr=kp(i);
        kd_curr=kd(j);
        % Quadrotor dynamics        
        dim=2;% spatial dimension (of positions and velocities)
        % Current implementation only supports dim=2 for quadrotors
        G_veh=define_G_quad_wrapped(dim,kp_curr,kd_curr);  
        [alpha_best(i,j),~]=bisection_exponent(G_veh,m,L,alpha_lims,cond_tol,cvx_tol,bisect_tol,multiplier_flag);        
    end
    if save_data==1
        save(save_path);
    end
end

