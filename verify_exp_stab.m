function [status,P]=verify_exp_stab(G,M,alpha,tol)
% This function runs the analysis LMI with cvx and return the status and
% the storage function matrix P. Inputs include the augmented plant G as
% Psi*[G;I], the fixed multiplier matrix M, exponent alpha and tolerance
% for the optimization

status=false;
[A,B,C,D]=ssdata(G);
n=size(A,1);
m=size(B,2);

cvx_begin sdp 
cvx_precision high
variable lambda
variable P(n,n) symmetric 

L1=[A'*P + P*A + 2*alpha*P,    P*B;
    B'*P,                      zeros(m)];
L2=lambda*[C';D']*M*[C,D];

LMI=L1+L2;

minimize 1; 
subject to:
P >= tol*eye(n);
LMI<=0;
lambda>=tol;
cvx_end 
if strcmp('Solved',cvx_status)
    status=true;
end
end