function [status,P]=verify_exp_stab_ZF(G,alpha,sec_1,sec_2,dim,tol)
% This function runs the analysis LMI with cvx and returns the status and
% the storage function matrix P. Inputs include the augmented plant G as
% Psi*[G;I], the fixed multiplier matrix M, exponent alpha and tolerance
% for the optimization
    status=false;
    [A_G,B_G,C_G,D_G]=ssdata(G);
    % Negative feedback(if gradient is the delta block, we want grad descent)
    B_G=-B_G;
    nx=size(A_G,1);
    [~,nu]=size(D_G);
    
    %Pi_ab=[sec_1,-1;-sec_2, 1];
    
    % Define the multiplier H
    a=-2;
    A_h=a;
    B_h=1;
    C_h=-a;
    D_h=0;
    
    % Define the tall psi basis. Stable causal filter
    
%     A_psi=A_h-2*alpha;
%     B_psi=1;
%     C_psi=-a;
    
    
%     % Define the augmented 
%     A_PSI=A_psi;
%     B_PSI=[0,B_psi]*Pi_ab;
%     C_PSI=[0;0;C_psi];
%     D_PSI=[eye(2);0,D_psi]*Pi_ab;

    % Define the augmented 
    A_PSI=A_h-2*alpha;
    B_PSI=B_h*[sec_2,-1];
    C_PSI=[-C_h;0;0;0];
    D_PSI=[ sec_2,-1;...
            -sec_1, 1;...
            sec_2,-1;...
            -sec_1, 1];
    
    PSI=ss(A_PSI,B_PSI,C_PSI,D_PSI);
    
    
%     g=1;c=1;
%     M=[0,g,-c;g,0,0;-c',0,0];

    M=[0,1;1,0];
    
    
    % Build Psi*[G;I] with Dynamic Multiplier
    PSI_GI=PSI*ss(A_G,B_G,[C_G;zeros(nu,nx)],[D_G;eye(nu)]);
    
    [A,B,C,D]=ssdata(PSI_GI);
    n=size(A,1);
    m=size(B,2);

    cvx_begin sdp 
    cvx_precision high
    variable lambda_1
    variable lambda_2
    variable P(n,n) symmetric 

    L1=[A'*P + P*A + 2*alpha*P,    P*B;
        B'*P,                      zeros(nu)];
    L2=[C';D']*[lambda_1*M,zeros(2);zeros(2),lambda_2*M]*[C,D];

    LMI=L1+L2;

    minimize 1; 
    subject to:
    P >= tol*eye(n);
    P <= 100*tol*eye(n);
    LMI<=-tol*eye(n+nu);
    lambda_1>=tol;
    lambda_2>=tol;
    cvx_end 
    if strcmp('Solved',cvx_status)
        status=true;
    end
end