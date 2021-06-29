function [status,X]=verify_exp_stab_FBCC(G,alpha,sec_1,sec_2,dim,tol)
% This function runs the analysis LMI with cvx and returns the status and
% the storage function matrix P. Inputs include the augmented plant G as
% Psi*[G;I], the fixed multiplier matrix M, exponent alpha and tolerance
% for the optimization
    status=false;
    [A,B,C,D]=ssdata(G);
    n=size(A,1);
    m=size(B,2);

    cvx_begin sdp 
    cvx_precision high    
    variable X(n,n) symmetric
    variable Q(dim,dim) symmetric
    variable S(dim,dim)
    variable R(dim,dim) symmetric

    % Multiplier class
    % Create the multiplier class
    d=[ sec_1,sec_1;
        sec_1,sec_2;
        sec_2,sec_1;
        sec_2,sec_2];
    E=@(delta)[eye(dim);diag(delta)];
    
    C1=E(d(1,:))'*[Q,S;S',R]*E(d(1,:));
    C2=E(d(2,:))'*[Q,S;S',R]*E(d(2,:));
    C3=E(d(3,:))'*[Q,S;S',R]*E(d(3,:));
    C4=E(d(4,:))'*[Q,S;S',R]*E(d(4,:));

    % FDI via KYP Lemma
    L1=[A'*X + X*A + 2*alpha*X,    X*B;
        B'*X,                      zeros(m)];
    L2=[C';D']*[Q,S;S',R]*[C,D];
    LMI=L1+L2;

    minimize 1; 
    subject to:
    R<= 0;
    C1>=0;
    C2>=0;
    C3>=0;
    C4>=0;
    
    X >= tol*eye(n);
    X <= 100*tol*eye(n);
    
    LMI<=0;    
    cvx_end 
    if strcmp('Solved',cvx_status)
        status=true;
    end
end