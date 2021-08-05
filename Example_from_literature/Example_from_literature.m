% This script is used to reproduce data for Examples from Literature
close all
clear
clc
addpath('..\analysis_scripts')
save_data=1;
save_path='.\';
%% Setup optimization

% Example 1 from Chen and Wen(not a vehicle): the optimal multiplier method
% for nonlinear robustness analysis
% G_veh=tf([3,3],[1,1,25,0,0]);

% Example 1 from Chen and Wen(not a vehicle): the optimal multiplier method
% for nonlinear robustness analysis
%G_veh=tf([1,1],[1,1,25,0,0]);

% Example from Scherer and Weiland
z1=0.4;
p1=-0.56984;
p2_r=-0.21508;
p2_i=1.3;
s=tf('s');
G_veh=-s*(s-z1)/((s-p1)*(s^2-2*(p2_r)*s+p2_r^2+p2_i^2));            

% Sector bounds
m=0; % lower bound on the sector
L=2:0.1:5;  % Upper bound on the sector
L=2:0.1:3;
n_L=length(L);

% Optimization properties
cvx_tol=1e-6;
bisect_tol=1e-2;
cond_tol=100000000;

% Multiplier class
% Select a multiplier class from the following choices
% 1. Circle criterion
% 2. Full block circle criterion
% 3. Zames Falb multipliers
multiplier_flag=6;
save_path=[save_path,'scherer_mult_flag_',num2str(multiplier_flag)];
alpha_lims=[0,10]; % Initial range for the bisection algorithm
alpha_best=zeros(1,n_L);
%% Run the bisection algorithm over a 2D grid (a X L)

for j=1:n_L
    if j>1 && alpha_best(1,j-1)==-1
        alpha_best(1,j)=-1;            
    else
        L_curr=L(j);            
        [alpha_best(1,j),~]=bisection_exponent(G_veh,m,L_curr,alpha_lims,cond_tol,cvx_tol,bisect_tol,multiplier_flag);
    end
end    

if save_data==1
    save(save_path);
end

