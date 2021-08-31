% This script is used to reproduce data for Example 7_20 from Scherer and
% Weiland's LMI notes. It is an exercise with an odd sector bounded
% non-linearity
close all
clear
clc
addpath('..\analysis_scripts')
save_data=0;
save_path='.\data_fourth_order_plant\';
%% Setup optimization
kp=-1;
n_p=length(kp);

% Sector bounds
m=1; % lower bound on the sector
L=1.1;  % Upper bound on the sector

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
save_path=[save_path,'mult_flag_',num2str(multiplier_flag)];
alpha_lims=[0,10]; % Initial range for the bisection algorithm
alpha_best=zeros(n_p,1);
%% Run the bisection algorithm over a 2D grid (a X L)
for i=1:n_p    
        kp_curr=kp(i);
        % Plant model
        G_veh=ss(kp_curr*tf([3,-3],[1,1,25,0,0]));      % Chen's example
        %G_veh=ss(kp_curr*tf([-1,1,0],[1,1,2,1]));      % Scherer's example
        
        [alpha_best(i,1),~]=bisection_exponent(G_veh,m,L,alpha_lims,cond_tol,cvx_tol,bisect_tol,multiplier_flag); 
end
if save_data==1
    save(save_path);
end

