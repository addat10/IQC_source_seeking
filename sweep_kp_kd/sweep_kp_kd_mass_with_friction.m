% This script is used to reproduce data for Example 7_20 from Scherer and
% Weiland's LMI notes. It is an exercise with an odd sector bounded
% non-linearity
close all
clear
clc
addpath('..\analysis_scripts')
save_data=1;
save_path='.\data_mass_with_friction\';
%% Setup optimization
kp=0.5:0.5:5;
kd=0.1:0.01:10;
n_p=length(kp);
n_d=length(kd);

% Sector bounds
m=1; % lower bound on the sector
L=20;  % Upper bound on the sector

% Optimization properties
cvx_tol=1e-3;
bisect_tol=1e-2;
cond_tol=100000000;

% Multiplier class
% Select a multiplier class from the following choices
% 1. Circle criterion
% 2. Full block circle criterion
% 3. Zames Falb multipliers
% 6. Zames Falb multipliers with basis
multiplier_flag=6;
save_path=[save_path,'mult_flag_non_causal_',num2str(multiplier_flag)];
alpha_lims=[0,10]; % Initial range for the bisection algorithm
alpha_best=zeros(n_p,n_d);
%% Run the bisection algorithm over a 2D grid (a X L)
for i=1:n_p
    for j=1:n_d        
        kp_curr=kp(i);
        kd_curr=kd(j);
        % Plant model
        G_veh=ss(tf([kp_curr],[1,kd_curr,0]));            
        [alpha_best(i,j),~]=bisection_exponent(G_veh,m,L,alpha_lims,cond_tol,cvx_tol,bisect_tol,multiplier_flag);        
    end    
end
if save_data==1
    save(save_path);
end

