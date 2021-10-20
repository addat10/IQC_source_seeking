% Investigate a lower bound on smallest eigen value of the grounded Laplacian
% M. Pirani and S. Sundaram, "On the Smallest Eigenvalue of Grounded 
% Laplacian Matrices," in IEEE Transactions on Automatic Control, vol. 61, 
% no. 2, pp. 509-514, Feb. 2016, doi: 10.1109/TAC.2015.2444191 
clear all
clc
close all
%% Generate a random Laplacian matrix
n=50; % total number of agents
m=10; % number of informed agents
sample=50;
link_prob=0.50;
topo='rand';
eig_min_vals=zeros(1,sample);
eig_lb=zeros(1,sample);
% Draw random Laplacians to test the bound
for i=1:sample
    A=gen_topology(n,link_prob,topo);
    D=A*ones(n,1);
    L=(diag(D)-A);
    q_vec=zeros(n,1);
    q_vec(1:m)=ones(m,1);
    Q=diag(q_vec);
    eig_L=eig(L);
    %% Compute the exact eigen value and the lower bound
    k=eig_L(2)/(4*sqrt(m));
    eig_vals=eig(L+k*Q);
    eig_min_vals(1,i)=eig_vals(1);
    eig_lb(1,i)=(k*m/n)*(1-2*k*sqrt(m)/eig_L(2));
    eig_ub(1,i)=(k*m/n);
end
%% Compare the lower bound with the actual smallest eigen value
figure()
plot(1:sample,eig_min_vals(1,:),'r')
hold on
plot(1:sample,eig_lb(1,:),'k')
plot(1:sample,eig_ub(1,:),'k')
