function [status,P]=verify_exp_stab_CC(G_veh,alpha,sec_1,sec_2,cond_tol,tol)
% This function runs the analysis LMI with cvx and return the status and
% the storage function matrix P. 
    
    [A_G,B_G,C_G,D_G]=ssdata(G_veh);
    
    % Negative feedback(if gradient is the delta block, we want grad descent)
    B_G=-B_G;
    
    nx=size(A_G,1);
    [ny,nu]=size(D_G); 
    
    % Define the static Multiplier
    A_PSI=0;
    B_PSI=zeros(1,ny+nu);
    C_PSI=zeros(ny+nu,1);
    D_PSI=[ sec_2*eye(ny),-1*eye(nu);...
            -sec_1*eye(ny), 1*eye(nu)];

    PSI=ss([],[],[],D_PSI);
    M=[ zeros(ny),eye(nu);...
        eye(ny),zeros(nu)];    
    
    % Build Psi*[G;I] with Dynamic Multiplier
    PSI_GI=PSI*ss(A_G,B_G,[C_G;zeros(nu,nx)],[D_G;eye(nu)]);
    
    % Run the exponenetial analysis LMI
    [status,P]=verify_exp_stab_alpha_CC_LMI(PSI_GI,M,alpha,cond_tol,tol);
    
end
function [status,P]=verify_exp_stab_alpha_CC_LMI(PSI_GI,M,alpha,cond_tol,tol)
 % This function runs the exp-stab analysis KYP Lemma LMI
    status=false;
    % Get state-space matrics
    [A,B,C,D]=ssdata(PSI_GI);    
    [n,nu]=size(B);

    cvx_begin sdp 
    cvx_precision high
    variable lambda
    variable P(n,n) symmetric 

    L1=[A'*P + P*A + 2*alpha*P,    P*B;
        B'*P,                      zeros(nu)];
    L2=[C';D']*lambda*M*[C,D];

    LMI=L1+L2;

    minimize 1; 
    subject to:
    P >= tol*eye(n);
    P <= cond_tol*tol*eye(n);
    LMI<=-tol*eye(n+nu);
    lambda>=tol;
    cvx_end 
    if strcmp('Solved',cvx_status)
        status=true;
    end
end