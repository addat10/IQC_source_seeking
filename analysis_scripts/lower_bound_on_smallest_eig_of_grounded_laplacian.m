% Investigate a lower bound on smallest eigen value of the grounded Laplacian
% M. Pirani and S. Sundaram, "On the Smallest Eigenvalue of Grounded 
% Laplacian Matrices," in IEEE Transactions on Automatic Control, vol. 61, 
% no. 2, pp. 509-514, Feb. 2016, doi: 10.1109/TAC.2015.2444191 
clear all
clc
close all
%% Generate a random Laplacian matrix
n=5;
link_prob=0.8;
topo='rand';
A=gen_topology(n,link_prob,topo);
D=A*ones(n,1);
L=(diag(D)-A);
Q=zeros(size(L));Q(1,1)=1;
eig_L=eig(L);
%% Compute the eigen values for grounded Laplacians with varying perturbations
k=0:0.1:2;
k_len=size(k,2);
eig_vals=zeros(n,k_len);
eig_lb=zeros(1,k_len);
for i=1:size(k,2)
    eig_vals(:,i)=eig(L+k(i)*Q);
    eig_lb(1,i)=(k(i)/n)*(1-2*k(i)/eig_L(2));
end
%% Compare the lower bound with the actual smallest eigen value
figure()
plot(k,eig_L(1)*ones(1,size(k,2)),'ro')
hold on
plot(k,eig_vals(1,:),'r')
% plot(eig_L(2)*ones(1,size(k,2)),'g-')
% plot(eig_vals(2,:),'g')
% plot(eig_L(3)*ones(1,size(k,2)),'b-')
% plot(eig_vals(3,:),'b')
plot(k,eig_lb,'k')
