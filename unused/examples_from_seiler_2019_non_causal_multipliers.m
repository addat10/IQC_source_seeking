close all
clear
clc
%% Vehicle Dynamics
% Example from Carasso, Seiler paper
a=-(0.4:0.1:0.4);
n_a=length(a);
%% Verify exponential stability: Analysis
addpath('.\analysis_scripts')
m=0;
L=0.1:0.1:0.1; 
n_L=length(L);
cvx_tol=1e-6;
bisect_tol=1e-2;
alpha_lims=[1e-6,10]; 
cond_tol=100000000;

% Multiplier class
% Select a multiplier class from the following choices
% 1. Circle criterion
% 2. Full block circle criterion
% 3. Zames Falb multipliers
alpha_best_CC=zeros(n_a,n_L);
alpha_best_ZFC=zeros(n_a,n_L);
alpha_best_ZFNC=zeros(n_a,n_L);

for i=1:n_a
    for j=1:n_L
        a_curr=a(i);
        L_curr=L(j);
        %G_veh=ss(tf([1,-a_curr,0],[1,1,2,1]));%%
        % G_veh=ss(tf([-1],[1,5]));
        G_veh=ss(tf([1,-2],[1,2,1,1]));
        % Circle criterion
        multiplier_flag=1;
        [alpha_best_CC(i,j),~]=bisection_exponent(G_veh,m,L_curr,alpha_lims,cond_tol,cvx_tol,bisect_tol,multiplier_flag);
        
        % ZF Causal
        multiplier_flag=3;
        [alpha_best_ZFC(i,j),~]=bisection_exponent(G_veh,m,L_curr,alpha_lims,cond_tol,cvx_tol,bisect_tol,multiplier_flag);
        
        % ZF Non-causal
        multiplier_flag=5;
        [alpha_best_ZFNC(i,j),~]=bisection_exponent(G_veh,m,L_curr,alpha_lims,cond_tol,cvx_tol,bisect_tol,multiplier_flag);
    end    
end
save('alpha_CC','alpha_best_CC','L');
save('alpha_ZFC','alpha_best_ZFC','L');
save('alpha_ZFNC','alpha_best_ZFNC','L');