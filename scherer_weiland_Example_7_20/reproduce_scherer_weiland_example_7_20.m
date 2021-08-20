% This script is used to reproduce data for Example 7_20 from Scherer and
% Weiland's LMI notes. It is an exercise with an odd sector bounded
% non-linearity
close all
clear
clc
addpath('..\analysis_scripts')
save_data=1;
save_path='.\data_non_odd_phi_alpha\';
%% Setup optimization
% Example from Scherer, Weiland LMI notes (not a vehicle)
a=0.4:0.2:1.8;
n_a=length(a);

% Sector bounds
m=0; % lower bound on the sector
L=0.1:0.01:5;  % Upper bound on the sector
n_L=length(L);

% Optimization properties
cvx_tol=1e-3;
bisect_tol=1e-2;
cond_tol=100000000;

% Multiplier class
% Select a multiplier class from the following choices
% 1. Circle criterion
% 2. Full block circle criterion
% 3. Zames Falb multipliers
multiplier_flag=6;
save_path=[save_path,'mult_flag_non_causal_',num2str(multiplier_flag)];
alpha_lims=[0,10]; % Initial range for the bisection algorithm
alpha_best=zeros(n_a,n_L);
%% Run the bisection algorithm over a 2D grid (a X L)
for i=1:n_a
    for j=1:n_L
        if j>1 && alpha_best(i,j-1)==-1
            alpha_best(i,j)=-1;            
        else
            a_curr=a(i);
            L_curr=L(j);
            % Plant model
            G_veh=ss(tf([1,-a_curr,0],[1,1,2,1]));            
            [alpha_best(i,j),~]=bisection_exponent(G_veh,m,L_curr,alpha_lims,cond_tol,cvx_tol,bisect_tol,multiplier_flag);
        end
    end    
end
if save_data==1
    save(save_path,'alpha_best','L','a');
end

