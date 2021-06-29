function [status,P]=verify_exp_stab_ZF(G,alpha,sec_1,sec_2,dim,tol)
% This function runs the analysis LMI with cvx and returns the status and
% the storage function matrix P. Inputs include the augmented plant G as
% Psi*[G;I], the fixed multiplier matrix M, exponent alpha and tolerance
% for the optimization
    status=false;
    [A_G,B_G,C_G,D_G]=ssdata(G);
    nx=size(A_G,1);
    [ny,nu]=size(D_G);
    
    Pi_ab=[sec_1,-1;-sec_2, 1];
    
    % Define the tall psi basis. Stable causal filter
    A_psi=-1;
    B_psi=1;
    C_psi=1;
    D_psi=0;
    
    % Define the augmented 
    A_PSI=A_psi;
    B_PSI=[0,B_psi]*Pi_ab;
    C_PSI=[0;0;C_psi];
    D_PSI=[eye(2);0,D_psi]*Pi_ab;
    
    PSI=ss(A_PSI,B_PSI,C_PSI,D_PSI);
    
    g=1;c=1;
    M=[0,g,-c;g,0,0;-c',0,0];
    
    
    % Build Psi*[G;I] with Dynamic Multiplier
    PSI_GI=PSI*ss(A_G,B_G,[C_G;zeros(nu,nx)],[D_G;eye(nu)]);
    
    [A,B,C,D]=ssdata(PSI_GI);
    n=size(A,1);
    m=size(B,2);

    cvx_begin sdp 
    cvx_precision high
    variable lambda
    variable P(n,n) symmetric 

    L1=[A'*P + P*A + 2*alpha*P,    P*B;
        B'*P,                      zeros(nu)];
    L2=lambda*[C';D']*M*[C,D];

    LMI=L1+L2;

    minimize 1; 
    subject to:
    P >= tol*eye(n);
    P <= 100*tol*eye(n);
    LMI<=0;
    lambda>=tol;
    cvx_end 
    if strcmp('Solved',cvx_status)
        status=true;
    end
end